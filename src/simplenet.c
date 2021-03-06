/*
  n is the number of values to return
  ncon is the number of connections each guy gets
  index[n][ncon] is the list of indices to which it connects
  weight[n][ncon] is the list of weights
  root is the address of lowest entry in the variable array

  To wit:
  for sparse type
  name(i) produces values[i]
  value[i] = sum(k=0..ncon-1)of(w[i][k]*root[index[i][k]])

ODE FILE CALL:

special f=sparse( n, ncon, w,index,rootname)
special f=sconv(type,n,ncon,w,rootname)
index,w are tables that are loaded by the "tabular"
command.
rootname here is the name of the root quantity.
type is either ep0
for conv.


conv0 conve convp type

value[i]=    sum(j=-ncon;j<=ncon;i++) {
			  k=i+j;
		0 type      if(k>=0 && k<n)
			  value[i]+=wgt[j+ncon]*rootname[k]
		e type k=abs(i+j); if(k<2n)
							   if(k>=n)k=n-k;
				   ...
		p type
		   k=mod(2*n+i+j,n)

for example  discretized diffusion
tabular dd % 3 -1 1 3*abs(t)-2
special diff=conv(even, 51, 2, dd,v0)
v[0..50]'=f(v[j],w[j])+d*diff([j])

another example
nnetwork

tabular wgt % 51 -25 25 .032*cos(.5*pi*t/25)
special stot=conv(0, 51, 25, wgt,v0)
v[0..50]'=f(v[j],w[j])-gsyn*stot([j])*(v([j])-vsyn)


last example -- random sparse network 51 cells 5 connections each

tabular w % 251 0 250 .2*rand(1)
tabular con % 251 0 250 flr(rand(1)*5)
special stot=sparse(51,5,w,con,v0)
v[0..50]'=f(v[j],w[j])-gsyn*stot([j])*(v([j])-vsyn)

more stuff:
special k= delmmult(n,m,w,tau,root)
		  w has n*m values and tau has n*m delays
		  n is the length of root (cols) and m is number of rows
		  k(i) = sum(j=0,n-1) w[i*n+j] delay(root(j),tau[i*n+j])
		  Note delays can only work on variables and not on
		  fixed  values
special k=delsparse(m,nc,w,index,tau,root)
		  m is the number of return values
		  nc is number of connections per node. Must be the sam
		  w is  n*nc is weights
		  index is  m*nc that has the indices of connections

		  tau is m*nc is list of delays

	   k[i]=sum( 0<=j < nc) w[j+nc*i]*delay(root[c[j+nc*i]],tau[j+nc*i])


special f=mmult(n,m,w,root)  -- note that here m is the number of return
						  values and n is the size of input root
f[j] = sum(i=0,n-1)of(w(i+n*j)*root(i)) j=0,..,m-1
special f=fmmult(n,m,w,root1,root2,fname)
f[j] = sum(i=0,n-1)of(w(i+n*j)*fname(root1[i],root2[j])
special f=fsparse( n, ncon,w,index,root1,root2,fname)
special f=fconv(type,n,ncon,w,root1,root2,fname)
sum(j=-ncon,ncon)w(j)*fname(root1[i+j],root2[i])
similarly for fsparse

special k=fftcon(type,n,wgt,v0)

uses the fft to convolve v0 with the wgt which must be of
length n if type=periodic
	   2n if type=0

special f=interp(type,n,root)
		this produces an interpolated value of x for
		x \in [0,n)
		type =0 for linear (all that works now)
		root0,....rootn-1  are values at the integers
		Like a table, but the integer values are variables

special f=findext(type,n,skip,root)
if type=1  mx(0)=maximum mx(1)=index
if type=-1 mx(2)=minimum mx(3)=index
if type=0 mx(0)=maximum mx(1)=index,mx(2)=minimum mx(3)=index


special ydot=import(soname,sofun,nret,root,w1,w2,...wm)

 run right hand side in C
 soname is shared object library say mlnet.so
 sofun  is shared object function
 nret is the number of return values
 root is the name of the first variable



sofun(int nret, int root, double *con, double *var, double *z[50],double *ydot)

*z[50] contains a list of pointers  z[0] -> w1, .... 50 is hard coded


NOTE that the user-defined parameters start at #6 and are in order
including derived parameters but XPP takes care of this so start at 0
*/
#include "simplenet.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "delay_handle.h"
#include "extra.h"
#include "fftn.h"
#include "form_ode.h"
#include "ggets.h"
#include "markov.h"
#include "parserslow.h"
#include "util/strutil.h"
#include "tabular.h"

/* --- Macros --- */
#define MAXW 50
#define EVEN 0
#define ZERO 1
#define PERIODIC 2
/* simple network stuff */

#define CONVE 0
#define CONV0 1
#define CONVP 2
#define SPARSE 3
#define FCONVE 10
#define FCONV0 11
#define FCONVP 12
#define FSPARSE 13
#define FFTCON0 4
#define FFTCONP 5
#define MMULT 6
#define FMMULT 7
#define GILLTYPE 25
#define INTERP 30
#define FINDEXT 35  /* find extrema in list of variables */
#define DEL_MUL 40  /* for delayed coupled networks  - global coupling */
#define DEL_SPAR 41 /* sparse with unequal in degree and delays  */
#define IMPORT  50 /* not really a network type   */


/* --- Types --- */
typedef struct {
	int type,ncon,n;
	char name[20];
	char soname[DEFAULT_STRING_LENGTH],sofun[DEFAULT_STRING_LENGTH];
	int root,root2;
	int f[20];
	int iwgt;
	int *gcom; /* for group commands */
	double *values,*weight,*index,*taud; /* for delays  */
	double *fftr,*ffti,*dr,*di;
	double *wgtlist[MAXW];
} NETWORK;

typedef struct {
	char name[20];
	int root,length,il,ir;
} VECTORIZER;


/* --- Forward declarations */
static void fft_conv(int it,int n,double *values,double *yy,double *fftr,double *ffti,double *dr,double *di);
static int g_namelist(char *s, char *root, int *flag, int *i1, int *i2);
static int getimpstr(char *in,int *i,char *out);
static int get_vector_info(char *str, char *name,int *root, int *length, int *il, int *ir);
static int gilparse(char *s, int *ind, int *nn);
static int import_error(void);
static void init_net(double *v, int n);
static int is_network(char *s);
static double net_interp(double x, int i);
static int parse_import(char *s,  char *soname,char *sofun,int *n, char *vname,int *m, char *tname[MAXW]);
static void update_fft(int ind);

/* extra.c */
void get_import_values(int n, double *ydot, char *soname, char *sofun,
					   int ivar, double *wgt[MAXW],
					   double *var, double *con);

/* --- Data --- */
static int n_vector=0;
static int n_network=0;
static VECTORIZER my_vec[MAXVEC];
static NETWORK my_net[MAXNET];


/* --- Functions --- */
int add_spec_fun(char *name, char *rhs) {
	int i,ind,elen;
	int type;
	int iwgt,itau,iind,ivar,ivar2;
	int ntype,ntot,ncon=0,ntab=0;
	char *str, *toksave;
	char junk[DEFAULT_STRING_LENGTH];
	char rootname[20],wgtname[20],tauname[20],indname[20];
	char root2name[20],fname[20];
	char sofun[DEFAULT_STRING_LENGTH],soname[DEFAULT_STRING_LENGTH],*tname[MAXW];
	type=is_network(rhs);
	if(type==0) {
		return 0;
	}
	plintf("type=%d \n",type);
	for(i=0;i<n_network;i++) {
		if(strcmp(name,my_net[i].name)==0){
			break;
		}
	}
	ind=i;
	if(ind>=n_network) {
		plintf(" No such name %s ?? \n",name);
		return 0;
	}
	switch(type) {
	case 1: /* convolution */
		strtok_r(rhs,"(", &toksave);
		str=strtok_r(NULL, ",", &toksave);
		ntype=-1;
		if(str[0]=='E') {
			ntype=CONVE;
		}
		if(str[0]=='0'||str[0]=='Z') {
			ntype=CONV0;
		}
		if(str[0]=='P') {
			ntype=CONVP;
		}
		if(ntype==-1) {
			plintf(" No such convolution type %s \n",str);
			return 0;
		}
		str=strtok_r(NULL, ",", &toksave);
		ntot=atoi(str);
		if(ntot<=0) {
			plintf(" %s must be positive int \n",str);
			return 0;
		}
		str=strtok_r(NULL, ",", &toksave);
		ncon=atoi(str);
		if(ncon<=0) {
			plintf(" %s must be positive int \n",str);
			return 0;
		}
		str=strtok_r(NULL, ",", &toksave);
		strcpy(wgtname,str);
		iwgt=find_lookup(wgtname);
		if(iwgt<0) {
			plintf("in network %s,  %s is not a table \n",
				   name,wgtname);
			return 0;
		}
		str=strtok_r(NULL, ")", &toksave);
		strcpy(rootname,str);
		ivar=get_var_index(rootname);
		if(ivar<0) {
			plintf(" In %s , %s is not valid variable\n",
				   name,rootname);
			return 0;
		}
		my_net[ind].values=(double *)malloc((ntot+1)*sizeof(double));
		init_net(my_net[ind].values,ntot);
		my_net[ind].weight=my_table[iwgt].y;
		my_net[ind].type=ntype;
		my_net[ind].root=ivar;
		my_net[ind].n=ntot;
		my_net[ind].ncon=ncon;
		plintf(" Added net %s type %d len=%d x %d using %s var[%d] \n",
			   name,ntype,ntot,ncon,wgtname,ivar);
		return 1;
		break;
	case 2: /* sparse */
		strtok_r(rhs,"(", &toksave);
		str=strtok_r(NULL, ",", &toksave);
		ntype=SPARSE;
		ntot=atoi(str);
		if(ntot<=0) {
			plintf(" %s must be positive int \n",str);
			return 0;
		}
		str=strtok_r(NULL, ",", &toksave);
		ncon=atoi(str);
		if(ncon<=0) {
			plintf(" %s must be positive int \n",str);
			return 0;
		}
		str=strtok_r(NULL, ",", &toksave);
		strcpy(wgtname,str);
		iwgt=find_lookup(wgtname);
		if(iwgt<0) {
			plintf("in network %s,  %s is not a table \n",
				   name,wgtname);
			return 0;
		}
		str=strtok_r(NULL, ",", &toksave);
		strcpy(indname,str);
		iind=find_lookup(indname);
		if(iind<0) {
			plintf("in network %s,  %s is not a table \n",
				   name,indname);
			return 0;
		}
		str=strtok_r(NULL, ")", &toksave);
		strcpy(rootname,str);
		ivar=get_var_index(rootname);
		if(ivar<0) {
			plintf(" In %s , %s is not valid variable\n",
				   name,rootname);
			return 0;
		}
		my_net[ind].values=(double *)malloc((ntot+1)*sizeof(double));
		init_net(my_net[ind].values,ntot);
		my_net[ind].weight=my_table[iwgt].y;
		my_net[ind].index=my_table[iind].y;

		my_net[ind].type=ntype;
		my_net[ind].root=ivar;
		my_net[ind].n=ntot;
		my_net[ind].ncon=ncon;
		plintf(" Added sparse %s len=%d x %d using %s var[%d]  and %s\n",
			   name,ntot,ncon,wgtname,ivar,indname );
		return 1;
		break;
	case 3: /* convolution */
		strtok_r(rhs,"(", &toksave);
		str=strtok_r(NULL, ",", &toksave);
		ntype=-1;
		if(str[0]=='E') {
			ntype=FCONVE;
		}
		if(str[0]=='0'||str[0]=='Z') {
			ntype=FCONV0;
		}
		if(str[0]=='P') {
			ntype=FCONVP;
		}
		if(ntype==-1) {
			plintf(" No such convolution type %s \n",str);
			return 0;
		}
		str=strtok_r(NULL, ",", &toksave);
		ntot=atoi(str);
		if(ntot<=0) {
			plintf(" %s must be positive int \n",str);
			return 0;
		}
		str=strtok_r(NULL, ",", &toksave);
		ncon=atoi(str);
		if(ncon<=0) {
			plintf(" %s must be positive int \n",str);
			return 0;
		}
		str=strtok_r(NULL, ",", &toksave);
		strcpy(wgtname,str);
		iwgt=find_lookup(wgtname);
		if(iwgt<0) {
			plintf("in network %s,  %s is not a table \n",
				   name,wgtname);
			return 0;
		}

		str=strtok_r(NULL, ",", &toksave);
		strcpy(rootname,str);
		ivar=get_var_index(rootname);
		if(ivar<0) {
			plintf(" In %s , %s is not valid variable\n",
				   name,rootname);
			return 0;
		}

		str=strtok_r(NULL, ",", &toksave);
		strcpy(root2name,str);
		ivar2=get_var_index(root2name);
		if(ivar2<0) {
			plintf(" In %s , %s is not valid variable\n",
				   name,root2name);
			return 0;
		}
		str=strtok_r(NULL, ")", &toksave);
		strcpy(fname,str);
		sprintf(junk,"%s(%s,%s)",fname,rootname,root2name);
		if(add_expr(junk,my_net[ind].f,&elen)) {
			plintf(" bad function %s \n",fname);
			return 0;
		}
		my_net[ind].values=(double *)malloc((ntot+1)*sizeof(double));
		init_net(my_net[ind].values,ntot);
		my_net[ind].weight=my_table[iwgt].y;
		my_net[ind].type=ntype;
		my_net[ind].root=my_net[ind].f[0]; /* this is strange - I am adding the compiled names */
		my_net[ind].root2=my_net[ind].f[1];
		my_net[ind].n=ntot;
		my_net[ind].ncon=ncon;
		plintf(" Added net %s type %d len=%d x %d using %s %s(var[%d],var[%d]) \n",
			   name,ntype,ntot,ncon,wgtname,fname,ivar,ivar2);
		return 1;
		break;
	case 4: /* sparse */
		strtok_r(rhs,"(", &toksave);
		str=strtok_r(NULL, ",", &toksave);
		ntype=FSPARSE;
		ntot=atoi(str);
		if(ntot<=0) {
			plintf(" %s must be positive int \n",str);
			return 0;
		}
		str=strtok_r(NULL, ",", &toksave);
		ncon=atoi(str);
		if(ncon<=0) {
			plintf(" %s must be positive int \n",str);
			return 0;
		}
		str=strtok_r(NULL, ",", &toksave);
		strcpy(wgtname,str);
		iwgt=find_lookup(wgtname);
		if(iwgt<0) {
			plintf("in network %s,  %s is not a table \n",
				   name,wgtname);
			return 0;
		}
		str=strtok_r(NULL, ",", &toksave);
		strcpy(indname,str);
		iind=find_lookup(indname);
		if(iind<0) {
			plintf("in network %s,  %s is not a table \n",
				   name,indname);
			return 0;
		}

		str=strtok_r(NULL, ",", &toksave);
		strcpy(rootname,str);
		ivar=get_var_index(rootname);

		if(ivar<0) {
			plintf(" In %s , %s is not valid variable\n",
				   name,rootname);
			return 0;
		}

		str=strtok_r(NULL, ",", &toksave);
		strcpy(root2name,str);
		ivar2=get_var_index(root2name);
		if(ivar2<0) {
			plintf(" In %s , %s is not valid variable\n",
				   name,root2name);
			return 0;
		}
		str=strtok_r(NULL, ")", &toksave);
		strcpy(fname,str);
		sprintf(junk,"%s(%s,%s)",fname,rootname,root2name);
		if(add_expr(junk,my_net[ind].f,&elen)) {
			plintf(" bad function %s \n",fname);
			return 0;
		}

		my_net[ind].values=(double *)malloc((ntot+1)*sizeof(double));
		init_net(my_net[ind].values,ntot);
		my_net[ind].weight=my_table[iwgt].y;
		my_net[ind].index=my_table[iind].y;

		my_net[ind].type=ntype;
		my_net[ind].root=my_net[ind].f[0]; /* this is strange - I am adding the compiled names */
		my_net[ind].root2=my_net[ind].f[1];
		my_net[ind].n=ntot;
		my_net[ind].ncon=ncon;
		plintf(" Sparse %s len=%d x %d using %s %s(var[%d],var[%d]) and %s\n",
			   name,ntot,ncon,wgtname,fname,ivar,ivar2,indname );
		return 1;
		break;
	case 5: /* fft convolution */
		strtok_r(rhs,"(", &toksave);
		str=strtok_r(NULL, ",", &toksave);
		ntype=-1;
		/* if(str[0]=='E')ntype=CONVE; */
		if(str[0]=='0'||str[0]=='Z') {
			ntype=FFTCON0;
		}
		if(str[0]=='P') {
			ntype=FFTCONP;
		}
		if(ntype==-1) {
			plintf(" No such fft convolution type %s \n",str);
			return 0;
		}
		str=strtok_r(NULL, ",", &toksave);
		ntot=atoi(str);
		if(ntot<=0) {
			plintf(" %s must be positive int \n",str);
			return 0;
		}
		str=strtok_r(NULL, ",", &toksave);
		strcpy(wgtname,str);
		iwgt=find_lookup(wgtname);
		if(iwgt<0) {
			plintf("in network %s,  %s is not a table \n",
				   name,wgtname);
			return 0;
		}
		ntab=get_lookup_len(iwgt);
		if(type==FFTCONP&&ntab<ntot) {
			plintf(" In %s, weight is length %d < %d \n",name,ntab,ntot);
			return 0;
		}
		if(type==FFTCON0&&ntab<(2*ntot)) {
			plintf(" In %s, weight is length %d < %d \n",name,ntab,2*ntot);
			return 0;
		}
		str=strtok_r(NULL, ")", &toksave);
		strcpy(rootname,str);
		ivar=get_var_index(rootname);
		if(ivar<0) {
			plintf(" In %s , %s is not valid variable\n",
				   name,rootname);
			return 0;
		}
		if(ntype==FFTCON0) {
			ncon=2*ntot;
		} else {
			ncon=ntot;
		}
		my_net[ind].fftr=(double *)malloc((ncon+2)*sizeof(double));
		my_net[ind].ffti=(double *)malloc((ncon+2)*sizeof(double));
		my_net[ind].dr=(double *)malloc((ncon+2)*sizeof(double));
		my_net[ind].di=(double *)malloc((ncon+2)*sizeof(double));
		my_net[ind].iwgt=iwgt;
		my_net[ind].values=(double *)malloc((ntot+1)*sizeof(double));
		init_net(my_net[ind].values,ntot);
		my_net[ind].weight=my_table[iwgt].y;
		my_net[ind].type=ntype;
		my_net[ind].root=ivar;
		my_net[ind].n=ntot;
		my_net[ind].ncon=ncon;
		update_fft(ind);
		plintf(" Added net %s type %d len=%d x %d using %s var[%d] \n",
			   name,ntype,ntot,ncon,wgtname,ivar);
		return 1;
		break;
	case 6:   /* MMULT    ntot=n,ncon=m  */
		strtok_r(rhs,"(", &toksave);
		str=strtok_r(NULL, ",", &toksave);
		ntype=MMULT;
		ntot=atoi(str);
		if(ntot<=0) {
			plintf(" %s must be positive int \n",str);
			return 0;
		}
		str=strtok_r(NULL, ",", &toksave);
		ncon=atoi(str);
		if(ncon<=0) {
			plintf(" %s must be positive int \n",str);
			return 0;
		}
		str=strtok_r(NULL, ",", &toksave);
		strcpy(wgtname,str);
		iwgt=find_lookup(wgtname);
		if(iwgt<0) {
			plintf("in network %s,  %s is not a table \n",
				   name,wgtname);
			return 0;
		}
		str=strtok_r(NULL, ")", &toksave);
		strcpy(rootname,str);
		ivar=get_var_index(rootname);
		if(ivar<0) {
			plintf(" In %s , %s is not valid variable\n",
				   name,rootname);
			return 0;
		}

		my_net[ind].values=(double *)malloc((ncon+1)*sizeof(double));
		init_net(my_net[ind].values,ncon);
		my_net[ind].weight=my_table[iwgt].y;

		my_net[ind].type=ntype;
		my_net[ind].root=ivar;
		my_net[ind].n=ncon;
		my_net[ind].ncon=ntot;
		plintf(" Added mmult %s len=%d x %d using %s var[%d]\n",
			   name,ntot,ncon,wgtname,ivar,indname );
		return 1;
		break;
	case 7:  /* FMMULT */
		strtok_r(rhs,"(", &toksave);
		str=strtok_r(NULL, ",", &toksave);
		ntype=FMMULT;
		ntot=atoi(str);
		if(ntot<=0) {
			plintf(" %s must be positive int \n",str);
			return 0;
		}
		str=strtok_r(NULL, ",", &toksave);
		ncon=atoi(str);
		if(ncon<=0) {
			plintf(" %s must be positive int \n",str);
			return 0;
		}
		str=strtok_r(NULL, ",", &toksave);
		strcpy(wgtname,str);
		iwgt=find_lookup(wgtname);
		if(iwgt<0) {
			plintf("in network %s,  %s is not a table \n",
				   name,wgtname);
			return 0;
		}
		str=strtok_r(NULL, ",", &toksave);
		strcpy(rootname,str);
		ivar=get_var_index(rootname);
		if(ivar<0) {
			plintf(" In %s , %s is not valid variable\n",
				   name,rootname);
			return 0;
		}
		str=strtok_r(NULL, ",", &toksave);
		strcpy(root2name,str);
		ivar2=get_var_index(root2name);
		if(ivar2<0) {
			plintf(" In %s , %s is not valid variable\n",
				   name,root2name);
			return 0;
		}
		str=strtok_r(NULL, ")", &toksave);
		strcpy(fname,str);
		sprintf(junk,"%s(%s,%s)",fname,rootname,root2name);
		if(add_expr(junk,my_net[ind].f,&elen)) {
			plintf(" bad function %s \n",fname);
			return 0;
		}
		my_net[ind].values=(double *)malloc((ncon+1)*sizeof(double));
		init_net(my_net[ind].values,ncon);
		my_net[ind].weight=my_table[iwgt].y;

		my_net[ind].type=ntype;
		my_net[ind].root=my_net[ind].f[0]; /* this is strange - I am adding the compiled names */
		my_net[ind].root2=my_net[ind].f[1];
		my_net[ind].n=ncon;
		my_net[ind].ncon=ntot;
		plintf(" Added fmmult %s len=%d x %d using %s %s(var[%d],var[%d])\n",
			   name,ntot,ncon,wgtname,fname,ivar,ivar2);
		return 1;
	case FINDEXT:
		strtok_r(rhs,"(", &toksave);
		str=strtok_r(NULL, ",", &toksave);
		ntype=atoi(str);
		if(ntype>1||ntype<(-1)) {
			plintf("In %s,  type =-1,0,1 not %s \n",
				   name,ntype);
			return 0;
		}
		str=strtok_r(NULL, ",", &toksave);
		ntot=atoi(str);
		if(ntot<=0) {
			plintf("In %s,  n>0 not %s \n",
				   name,ntot);
			return 0;
		}
		str=strtok_r(NULL, ",", &toksave);
		ncon=atoi(str);
		if(ncon<=0) {
			plintf("In %s,  skip>=1 not %s \n",
				   name,ncon);
			return 0;
		}
		str=strtok_r(NULL, ")", &toksave);
		strcpy(rootname,str);
		ivar=get_var_index(rootname);
		if(ivar<0) {
			plintf(" In %s , %s is not valid variable\n",
				   name,rootname);
			return 0;
		}
		my_net[ind].values=(double *)malloc(6*sizeof(double));
		my_net[ind].type=FINDEXT;
		my_net[ind].root=ivar;
		my_net[ind].n=ntot;
		my_net[ind].ncon=ncon;
		my_net[ind].iwgt=ntype;
		plintf(" Added findextr %s: type=%d len=%d  skip= %d using var[%d] \n",
			   name,ntype,ntot,ncon,ivar);
		return 1;
	case 30:
		strtok_r(rhs,"(", &toksave);
		str=strtok_r(NULL, ",", &toksave);
		ivar=atoi(str);
		my_net[ind].type=INTERP;
		my_net[ind].iwgt=ivar;
		str=strtok_r(NULL, ",", &toksave);
		ivar=atoi(str);
		if(ivar<1) {
			plintf("Need more than 1 entry for interpolate\n");
			return 0;
		}
		my_net[ind].n=ivar; /* # entries in array */
		str=strtok_r(NULL, ")", &toksave);
		strcpy(rootname,str);
		ivar=get_var_index(rootname);
		if(ivar<0) {
			plintf(" In %s , %s is not valid variable\n",
				   name,rootname);
			return 0;
		}
		my_net[ind].root=ivar;
		plintf("Added interpolator %s length %d on %s \n",name,my_net[ind].n,rootname);
		return 1;
	case IMPORT:
		ntype=IMPORT;
		for(i=0;i<MAXW;i++) {
			tname[i]=(char *)malloc(25);
		}
		parse_import(rhs,soname,sofun,&ncon,rootname,&ntab,tname);
		my_net[ind].values=(double *)malloc((ncon+1)*sizeof(double));
		init_net(my_net[ind].values,ncon);
		my_net[ind].n=ncon;
		ivar=get_var_index(rootname);
		if(ivar<0) {
			plintf(" In %s , %s is not valid variable\n",
				   name,rootname);
			return 0;
		}
		strcpy(my_net[ind].soname,soname);
		strcpy(my_net[ind].sofun,sofun);
		my_net[ind].root=ivar;
		my_net[ind].type=ntype;
		my_net[ind].ncon=0;
		for(i=0;i<ntab;i++) {
			iwgt=find_lookup(tname[i]);
			plintf("Found %s\n",tname[i]);
			if(iwgt<0) {
				plintf("in network %s,  %s is not a table \n",
					   name,wgtname);
				return 0;
			}
			my_net[ind].wgtlist[i]=my_table[iwgt].y;
		}
		for(i=0;i<MAXW;i++) {
			free(tname[i]);
		}
		plintf(" Added import %s len=%d  with %s %s var[%d] %d weights\n",
			   name,my_net[ind].n,soname,sofun,ivar,ntab );
		return 1;
	case DEL_MUL:
		strtok_r(rhs,"(", &toksave);
		str=strtok_r(NULL, ",", &toksave);
		ntype=DEL_MUL;
		ntot=atoi(str);
		if(ntot<=0) {
			plintf(" %s must be positive int \n",str);
			return 0;
		}
		str=strtok_r(NULL, ",", &toksave);
		ncon=atoi(str);
		if(ncon<=0) {
			plintf(" %s must be positive int \n",str);
			return 0;
		}
		str=strtok_r(NULL, ",", &toksave);
		strcpy(wgtname,str);
		iwgt=find_lookup(wgtname);
		if(iwgt<0) {
			plintf("in network %s,  %s is not a table \n",
				   name,wgtname);
			return 0;
		}
		str=strtok_r(NULL, ",", &toksave);
		strcpy(tauname,str);
		itau=find_lookup(tauname);
		if(itau<0) {
			plintf("in network %s,  %s is not a table \n",
				   name,tauname);
			return 0;
		}
		str=strtok_r(NULL, ")", &toksave);
		strcpy(rootname,str);
		ivar=get_var_index(rootname);
		if(ivar<0) {
			plintf(" In %s , %s is not valid variable\n",
				   name,rootname);
			return 0;
		}
		my_net[ind].values=(double *)malloc((ncon+1)*sizeof(double));
		init_net(my_net[ind].values,ncon);
		my_net[ind].weight=my_table[iwgt].y;
		my_net[ind].taud=my_table[itau].y;

		my_net[ind].type=ntype;
		my_net[ind].root=ivar;
		my_net[ind].n=ncon;
		my_net[ind].ncon=ntot;
		plintf(" Added del_mul %s len=%d x %d using %s var[%d] with delay %s\n",
			   name,ntot,ncon,wgtname,ivar,indname,tauname );
		NDELAYS=1;
		return 1;
		break;
		return 0;
	case DEL_SPAR:
		strtok_r(rhs,"(", &toksave);
		str=strtok_r(NULL, ",", &toksave);
		ntype=DEL_SPAR;
		ntot=atoi(str);
		if(ntot<=0) {
			plintf(" %s must be positive int \n",str);
			return 0;
		}
		str=strtok_r(NULL, ",", &toksave);
		ncon=atoi(str);
		if(ncon<=0) {
			plintf(" %s must be positive int \n",str);
			return 0;
		}
		str=strtok_r(NULL, ",", &toksave);
		strcpy(wgtname,str);
		iwgt=find_lookup(wgtname);
		if(iwgt<0) {
			plintf("in network %s,  %s is not a table \n",
				   name,wgtname);
			return 0;
		}
		str=strtok_r(NULL, ",", &toksave);
		strcpy(indname,str);
		iind=find_lookup(indname);
		if(iind<0) {
			plintf("in network %s,  %s is not a table \n",
				   name,indname);
			return 0;
		}
		str=strtok_r(NULL, ",", &toksave);
		strcpy(tauname,str);
		itau=find_lookup(tauname);
		if(itau<0) {
			plintf("in network %s,  %s is not a table \n",
				   name,tauname);
			return 0;
		}
		str=strtok_r(NULL, ")", &toksave);
		strcpy(rootname,str);
		ivar=get_var_index(rootname);
		if(ivar<0) {
			plintf(" In %s , %s is not valid variable\n",
				   name,rootname);
			return 0;
		}
		my_net[ind].values=(double *)malloc((ntot+1)*sizeof(double));
		init_net(my_net[ind].values,ntot);
		my_net[ind].weight=my_table[iwgt].y;
		my_net[ind].index=my_table[iind].y;
		my_net[ind].taud=my_table[itau].y;
		my_net[ind].type=ntype;
		my_net[ind].root=ivar;
		my_net[ind].n=ntot;
		my_net[ind].ncon=ncon;
		plintf(" Added sparse %s len=%d x %d using %s var[%d]  and %s with dely %s\n",
			   name,ntot,ncon,wgtname,ivar,indname,tauname );
		NDELAYS=1;
		return 1;
		break;
		return 0;
	case 10:
		/*
	   z=GILL(meth,rxn list)
	   e.g
	   z=GILL(meth,r{1-15})
	   GILL is different -
	   iwgt=evaluation method - 0 is standard
								1 - tau-leap
	   root=number of reactions
	   values[0]=time of next reaction
	   values[1..root]=number of times this rxn took place
	   gcom contains list of all the fixed holding the reactions
	*/
		strtok_r(rhs,"(", &toksave);
		str=strtok_r(NULL, ",", &toksave);
		ivar=atoi(str);
		str=strtok_r(NULL, ")", &toksave);
		my_net[ind].type=GILLTYPE;
		if(ivar>0) {
			plintf(" Tau leaping not implemented yet. Changing to 0\n");
			ivar=0;
		}
		my_net[ind].iwgt=ivar;
		my_net[ind].gcom=(int *)malloc(1000*sizeof(int));
		if(gilparse(str,my_net[ind].gcom,&ivar2)==0) {
			return 0;
		}
		my_net[ind].root=ivar2;
		my_net[ind].n=ivar2+1;
		my_net[ind].ncon=-1;
		my_net[ind].values=(double *)malloc((ivar2+2)*sizeof(double));
		plintf("Added gillespie chain with %d reactions \n",ivar2);
		return 1;
	}
	return 0;
}


void add_special_name(char *name, char *rhs) {
	if(is_network(rhs)) {
		plintf(" netrhs = |%s| \n",rhs);
		if(n_network>=MAXNET) {
			return;
		}
		strcpy(my_net[n_network].name,name);
		add_net_name(n_network,name);
		n_network++;
	} else {
		plintf(" No such special type ...\n");
	}
}


int add_vectorizer(char *name,char *rhs) {
	int i,ivar,il,ir;
	int ind;
	int len;
	int flag;

	for(i=0;i<n_vector;i++) {
		if(strcmp(name,my_vec[i].name)==0) {
			break;
		}
	}
	ind=i;
	flag=get_vector_info(rhs,name,&ivar,&len,&il,&ir);

	if(flag==0) {
		return 0;
	}
	my_vec[ind].root=ivar;
	my_vec[ind].length=len;
	my_vec[ind].il=il;
	my_vec[ind].ir=ir;
	plintf("adding vector %s based on variable %d of length %d ends %d %d\n",
		   name,ivar,len,il,ir);
	return 1;
}


void add_vectorizer_name(char *name, char *rhs) {
	if(n_vector>=MAXVEC) {
		plintf("Too many vectors \n");
		exit(0);
	}
	strcpy(my_vec[n_vector].name,name);
	if(add_vector_name( n_vector,name)) {
		exit(0);
	}
	n_vector++;

}


void eval_all_nets(void) {
	int i;
	for(i=0;i<n_network;i++) {
		evaluate_network(i);
	}
}


void evaluate_network(int ind) {
	int i,j,k,ij;
	int imin,imax;
	double ymin,ymax;
	int skip;
	int mmt;
	int in0;
	double sum,z;
	int n=my_net[ind].n,*f;
	int ncon=my_net[ind].ncon;
	double *w,*y,*cc,*values,*tau;
	int twon=2*n,root=my_net[ind].root,root2=my_net[ind].root2;
	cc=my_net[ind].index;
	w=my_net[ind].weight;
	values=my_net[ind].values;
	/*  y=&variables[my_net[ind].root]; */
	switch(my_net[ind].type) {
	case FINDEXT:
		y=&variables[root];
		mmt=my_net[ind].iwgt;
		skip=ncon;
		imax=0;
		ymax=y[0];
		imin=0;
		ymin=y[0];
		i=0;

		if(mmt==1||mmt==0) { /* get max  */
			while(i<n) {
				if(y[i]>ymax) {
					imax=i;
					ymax=y[i];
				}
				i+=skip;
			}
		}
		if(mmt==(-1)||mmt==0) { /* get min */
			i=0;
			while(i<n) {
				if(y[i]<ymin) {
					imin=i;
					ymin=y[i];
				}
				i+=skip;
			}
		}
		values[1]=(double)imax;
		values[0]=ymax;
		values[3]=(double)imin;
		values[2]=ymin;
		break;
	case INTERP: /* do nothing! */
		break;
	case GILLTYPE:
		if(my_net[ind].ncon==-1&&my_net[ind].iwgt>0) {
			my_net[ind].weight=(double *)malloc(my_net[ind].root*NODE*sizeof(double));
			make_gill_nu(my_net[ind].weight,NODE,my_net[ind].root,my_net[ind].values);
			my_net[ind].ncon=0;
		}
		one_gill_step(my_net[ind].iwgt,my_net[ind].root,my_net[ind].gcom,my_net[ind].values);
		break;
	case CONVE:
		y=&variables[root];
		for(i=0;i<n;i++) {
			sum=0.0;
			for(j=-ncon;j<=ncon;j++) {
				k=abs(i+j);
				if(k<twon) {
					if(k>=n) {
						k=abs(twon-2-k);
					}
					sum+=(w[j+ncon]*y[k]);
				}
			}
			values[i]=sum;
		}
		break;
	case CONV0:
		y=&variables[root];
		for(i=0;i<n;i++) {
			sum=0.0;
			for(j=-ncon;j<=ncon;j++) {
				k=i+j;
				if(k<n&&k>=0) {
					sum+=(w[j+ncon]*y[k]);
				}
			}
			values[i]=sum;
		}
		break;
	case CONVP:
		y=&variables[root];
		for(i=0;i<n;i++) {
			sum=0.0;
			for(j=-ncon;j<=ncon;j++) {
				k=((twon+i+j)%n);
				sum+=(w[j+ncon]*y[k]);
			}
			values[i]=sum;
		}
		break;
	case FFTCONP:
		y=&variables[root];
		fft_conv(0,n,values,y,my_net[ind].fftr,my_net[ind].ffti,my_net[ind].dr,my_net[ind].di);
		break;
	case FFTCON0:
		y=&variables[root];
		fft_conv(1,n,values,y,my_net[ind].fftr,my_net[ind].ffti,my_net[ind].dr,my_net[ind].di);
		break;
	case IMPORT:
		get_import_values(n,values,my_net[ind].soname,my_net[ind].sofun,my_net[ind].root,my_net[ind].wgtlist,variables,&constants[6]);
		break;
	case DEL_MUL:
		tau=my_net[ind].taud;
		in0=my_net[ind].root;
		for(j=0;j<n;j++) {
			sum=0.0;
			for(i=0;i<ncon;i++) {
				ij=j*ncon+i;
				sum+=(w[ij]*get_delay(i+in0,tau[ij]));
			}
			values[j]=sum;
		}
		break;
	case MMULT:
		y=&variables[root];
		for(j=0;j<n;j++) {
			sum=0.0;
			for(i=0;i<ncon;i++) {
				ij=j*ncon+i;
				sum+=(w[ij]*y[i]);
			}
			values[j]=sum;
		}
		break;
	case DEL_SPAR:
		tau=my_net[ind].taud;
		in0=my_net[ind].root;
		for(i=0;i<n;i++) {
			sum=0.0;
			for(j=0;j<ncon;j++) {
				ij=i*ncon+j;
				k=(int)cc[ij];
				if(k>=0) {
					sum+=(w[ij]*get_delay(k+in0,tau[ij]));
				}
			}
			values[i]=sum;
		}
		break;
	case SPARSE:
		y=&variables[root];
		for(i=0;i<n;i++) {
			sum=0.0;
			for(j=0;j<ncon;j++) {
				ij=i*ncon+j;
				k=(int)cc[ij];
				if(k>=0) {
					sum+=(w[ij]*y[k]);
				}
			}
			values[i]=sum;
		}
		break;
		/*     f stuff  */
	case FCONVE:
		f=my_net[ind].f;
		for(i=0;i<n;i++) {
			sum=0.0;
			f[1]=root+i;
			for(j=-ncon;j<=ncon;j++) {
				k=abs(i+j);
				if(k<twon) {
					if(k>=n) {
						k=abs(twon-2-k);
					}
					f[0]=root2+k;
					z=evaluate(f);
					sum+=(w[j+ncon]*z);
				}
			}
			values[i]=sum;
		}
		break;
	case FCONV0:
		f=my_net[ind].f;
		for(i=0;i<n;i++) {
			sum=0.0;
			f[1]=root+i;
			for(j=-ncon;j<=ncon;j++) {
				k=i+j;
				if(k<n&&k>=0) {
					f[0]=root2+k;
					z=evaluate(f);
					sum+=(w[j+ncon]*z);
				}
			}
			values[i]=sum;
		}
		break;
	case FCONVP:
		f=my_net[ind].f;
		for(i=0;i<n;i++) {
			f[1]=root+i;
			sum=0.0;
			for(j=-ncon;j<=ncon;j++) {
				k=((twon+i+j)%n);
				f[0]=root2+k;
				z=evaluate(f);
				sum+=(w[j+ncon]*z);
			}
			values[i]=sum;
		}
		break;
	case FSPARSE:
		f=my_net[ind].f;
		for(i=0;i<n;i++) {
			f[1]=root+i;
			sum=0.0;
			for(j=0;j<ncon;j++) {
				ij=i*ncon+j;
				k=(int)cc[ij];
				if(k>=0) {
					f[0]=root2+k;
					z=evaluate(f);
					sum+=(w[ij]*z);
				}
			}
			values[i]=sum;
		}
		break;
	case FMMULT:
		f=my_net[ind].f;
		for(j=0;j<n;j++) {
			f[1]=root+j;
			sum=0.0;
			for(i=0;i<ncon;i++) {
				ij=j*ncon+i;
				f[0]=root2+i;
				z=evaluate(f);
				sum+=(w[ij]*z);
			}
			values[j]=sum;
		}
		break;
	}
}


double network_value(double x, int i) {
	int j=(int)x;
	if(my_net[i].type==INTERP) {
		return net_interp(x,i);
	}
	if(j>=0&&j<my_net[i].n) {
		return my_net[i].values[j];
	}
	return 0.0;
}


void update_all_ffts(void) {
	int i;

	for(i=0;i<n_network;i++) {
		if(my_net[i].type==FFTCON0||my_net[i].type==FFTCONP) {
			update_fft(i);
		}
	}
}


double vector_value(double x, int i) {
	int il=my_vec[i].il,ir=my_vec[i].ir,n=my_vec[i].length,k=(int)x;
	int root=my_vec[i].root;
	if(il==PERIODIC) {
		return variables[root+((k+n)%n)];
	}
	if(k<0) {
		if(il==ZERO) {
			return 0.0;
		}
		return variables[root-k-1];
	} else if((k>=0)&&(k<n))  {
		return variables[root+k];
	} else {
		if(ir==ZERO) {
			return 0.0;
		}
		return variables[2*n-k-1+root];
	}
}


/* --- Static functions --- */
static void fft_conv(int it,int n,double *values,double *yy,double *fftr,double *ffti,double *dr,double *di) {
	int i;
	int dims[2];
	double x,y;
	int n2=2*n;
	switch(it) {
	case 0:
		dims[0]=n;
		for(i=0;i<n;i++) {
			di[i]=0.0;
			dr[i]=yy[i];
		}
		fftn(1,dims,dr,di,1,-2.0);
		for(i=0;i<n;i++) {
			x=dr[i]*fftr[i]-di[i]*ffti[i];
			y=dr[i]*ffti[i]+di[i]*fftr[i];
			dr[i]=x;
			di[i]=y;
		}

		fftn(1,dims,dr,di,-1,-2.0);
		for(i=0;i<n;i++) {
			values[i]=dr[i];
		}
		return;
	case 1:
		dims[0]=n2;
		for(i=0;i<n2;i++) {
			di[i]=0.0;
			if(i<n) {
				dr[i]=yy[i];
			} else {
				dr[i]=0.0;
			}
		}
		fftn(1,dims,dr,di,1,-2.0);
		for(i=0;i<n2;i++) {
			x=dr[i]*fftr[i]-di[i]*ffti[i];
			y=dr[i]*ffti[i]+di[i]*fftr[i];
			dr[i]=x;
			di[i]=y;
		}
		fftn(1,dims,dr,di,-1,-2.0);
		for(i=0;i<n;i++) {
			values[i]=dr[i];
		}
		return;
	}
}


/* plucks info out of  xxx{aa-bb}  or returns string */
static int g_namelist(char *s,char *root,int *flag,int *i1,int*i2) {
	int i,n=strlen(s),ir=-1,j=0;
	char c,num[20];
	*flag=0;
	for(i=0;i<n;i++) {
		if(s[i]=='{') {
			ir=i;
		}
	}
	if(ir<0) {
		strcpy(root,s);
		return 1;
	}
	for(i=0;i<ir;i++) {
		root[i]=s[i];
	}
	root[ir]=0;
	*flag=1;
	j=0;
	for(i=ir+1;i<n;i++) {
		c=s[i];
		if(c=='-') {
			break;
		}
		num[j]=c;
		j++;
	}
	if(i==n) {
		plintf("Illegal syntax %s\n",s);
		return 0;
	}
	num[j]=0;
	*i1=atoi(num);
	ir=i+1;
	j=0;
	for(i=ir;i<n;i++) {
		c=s[i];
		if(c=='}') {
			break;
		}
		num[j]=c;
		j++;
	}
	num[j]=0;
	*i2=atoi(num);
	return 1;
}


static int getimpstr(char *in,int *i,char *out) {
	int j=0;
	int done=1;
	char c;
	int k=0;
	while(done>0)	{
		c=in[*i];
		if((c==',')||(c==')')) {
			out[j]=0;
			*i=*i+1;
			done=0;
			if(c==')') {
				k=1;
			} else {
				k=0;
			}
		} else {
			out[j]=c;
			j++;
			*i=*i+1;
		}
	}
	return k;
}


static int get_vector_info(char *str, char *name,int *root, int *length, int *il, int *ir) {
	int i=0;
	int ivar;
	int n=strlen(str);
	char c;
	int j;
	char temp[100];

	de_space(str);
	for(i=0;i<n;i++) {
		if(str[i]=='(') {
			break;
		}
	}
	i++;
	j=0;
	while(1) {
		c=str[i];
		if(c==',') {
			i++;
			break;
		}
		temp[j]=c;
		i++;
		j++;
	}
	temp[j]=0;
	ivar=get_var_index(temp);
	if(ivar<0) {
		plintf(" In vector %s , %s is not valid variable\n",
			   name,temp);
		return 0;
	}
	*root=ivar;
	j=0;
	while(1) {
		c=str[i];
		if(c==',') {
			i++;
			break;
		}
		temp[j]=c;
		i++;
		j++;
	}
	temp[j]=0;
	n=atoi(temp);

	*length=n;
	*il=PERIODIC;
	if(str[i]=='e'|| str[i]=='E') {
		*il=EVEN;
	}
	if(str[i]=='z'|| str[i]=='Z') {
		*il=ZERO;
	}
	i++;
	i++;
	*ir=PERIODIC;
	if(str[i]=='e'|| str[i]=='E') {
		*ir=EVEN;
	}
	if(str[i]=='z'|| str[i]=='Z') {
		*ir=ZERO;
	}
	return 1;
}


/* parsing stuff to get gillespie code quickly */
static int gilparse(char *s,int *ind,int *nn) {
	int i=0,n=strlen(s);
	char piece[50],b[20],bn[25],c;
	int i1,i2,jp=0,f;
	int k=0,iv;
	int id,m;
	plintf("s=|%s|",s);
	while(1) {
		c=s[i];
		if(c==','||i>(n-1)) {
			piece[jp]=0;
			if(g_namelist(piece,b,&f,&i1,&i2)==0) {
				printf("Bad gillespie list %s\n",s);
				return 0;
			}
			if(f==0) {
				plintf("added %s\n",b);
				iv=get_var_index(b);
				if(iv<0) {
					plintf("No such name %s\n",b);
					return 0;
				}
				ind[k]=iv;
				k++;
			} else {
				plintf("added %s{%d-%d}\n",b,i1,i2);
				m=i2-i1+1;
				for(id=0;id<m;id++) {
					sprintf(bn,"%s%d",b,id+i1);
					iv=get_var_index(bn);
					if(iv<0) {
						plintf("No such name %s\n",bn);
						return 0;
					}
					ind[k]=iv;
					k++;
				}
			}
			if(i>(n-1)) {
				*nn=k;
				return 1;
			}
			jp=0;
		} else {
			piece[jp]=c;
			jp++;
		}
		i++;
	}
	*nn=k;
	return 1;
}


static int import_error(void) {
	printf("k=import(soname,sofun,nret,var0,w1,...,wm)");
	return 0;
}


static void init_net(double *v,int n){
	int i;
	for(i=0;i<n;i++) {
		v[i]=0.0;
	}
}


static int is_network(char *s) {
	/*int n;
  */
	de_space(s);
	strupr(s);
	/* n=strlen(s); Not used*/
	if(s[0]=='C' &&s[1]=='O' &&s[2]=='N' && s[3]=='V') {
		return 1;
	}
	if(s[0]=='S' &&s[1]=='P' &&s[2]=='A' && s[3]=='R') {
		return 2;
	}
	if(s[0]=='F'&&s[1]=='C' &&s[2]=='O' &&s[3]=='N' && s[4]=='V'){
		return 3;
	}
	if(s[0]=='F' && s[1]=='S' &&s[2]=='P' &&s[3]=='A' && s[4]=='R') {
		return 4;
	}
	if(s[0]=='F' && s[1]=='F' && s[2]=='T' && s[3]=='C' ) {
		return 5;
	}
	if(s[0]=='M' &&s[1]=='M' &&s[2]=='U' && s[3]=='L') {
		return 6;
	}
	if(s[0]=='F'&& s[1]=='M' &&s[2]=='M' &&s[3]=='U' && s[4]=='L') {
		return 7;
	}
	if(s[0]=='G'&& s[1]=='I' &&s[2]=='L' &&s[3]=='L') {
		return 10;
	}
	if(s[0]=='I' && s[1]=='N' &&s[2]=='T' && s[3]=='E' && s[4]=='R') {
		return INTERP;
	}
	if(s[0]=='F' && s[1]=='I' &&s[2]=='N' && s[3]=='D' && s[4]=='E') {
		return FINDEXT;
	}
	if(s[0]=='D' &&s[1]=='E' &&s[2]=='L' && s[3]=='M') {
		return DEL_MUL;
	}
	if(s[0]=='D' &&s[1]=='E' &&s[2]=='L' && s[3]=='S') {
		return DEL_SPAR;
	}
	if(s[0]=='I' &&s[1]=='M' &&s[2]=='P' && s[3]=='O') {
		return IMPORT;
	}
	return 0;
}


static double net_interp(double x, int i) {
	int jlo=(int)x;
	double *y;
	int n=my_net[i].n;
	double dx=x-(double)jlo;
	y=&variables[my_net[i].root];
	if(jlo<0 || jlo>(n-1)) {
		return 0.0; /* out of range */
	}
	return (1-dx)*y[jlo]+dx*y[jlo+1];
}


static int parse_import(char *s,  char *soname,char *sofun,int *n, char *vname,int *m, char *tname[MAXW]) {
	char temp[DEFAULT_STRING_LENGTH];
	int j;
	char c;
	int i=0;
	int done=1;

	while(done>0) {
		c=s[i];
		i++;
		if(c=='(') {
			done=0;
		}
	}

	j=getimpstr(s,&i,temp);
	strcpy(soname,temp);
	if(j==1) {
		return(import_error());
	}
	j=getimpstr(s,&i,temp);
	strcpy(sofun,temp);
	if(j==1) {
		return(import_error());
	}
	j=getimpstr(s,&i,temp);
	*n=atoi(temp);
	if(j==1||*n<=0) {
		return(import_error());
	}
	j=getimpstr(s,&i,temp);
	strcpy(vname,temp);
	/*  plintf("%s %s %d %s\n",soname,sofun,*n,vname); */
	*m=0;
	if(j==1) {
		printf("No weights....\n");
		return(1);
	}
	done=1;
	while(done>0) {
		j=getimpstr(s,&i,temp);
		strcpy(tname[*m],temp);
		*m=*m+1;
		if(j==1) {
			done=0;
		}
	}
	return 1;
}


/*
 tabular weights are of size 2k+1
 and go from -k ... k
 for FFT's
 they are reordered as follows
 fftr[i]=wgt[i+k] i = 0.. k
 fftr[i+k]=wgt[i] i=1 .. k-1
*/
static void update_fft(int ind) {
	int i;
	int dims[2];
	double *w=my_net[ind].weight;
	double *fftr=my_net[ind].fftr;
	double *ffti=my_net[ind].ffti;
	int n,n2;
	int type=my_net[ind].type;
	if(type==FFTCONP) {
		n=my_net[ind].n;
		n2=n/2;
		for(i=0;i<n;i++) {
			ffti[i]=0.0;
		}
		for(i=0;i<=n2;i++) {
			fftr[i]=w[i+n2];
		}
		for(i=0;i<n2;i++) {
			fftr[n2+i+1]=w[i];
		}
		dims[0]=n;
		fftn(1,dims,fftr,ffti,1,1.);
	}
	if(type==FFTCON0) {
		n=2*my_net[ind].n;
		n2=n/2;
		for(i=0;i<n;i++) {
			ffti[i]=0.0;
		}
		for(i=0;i<=n2;i++) {
			fftr[i]=w[i+n2];
		}
		for(i=1;i<n2;i++) {
			fftr[n2+i]=w[i];
		}
		dims[0]=n;
		fftn(1,dims,fftr,ffti,1,1.);
	}
}
