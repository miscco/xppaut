#include "gear.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "abort.h"
#include "eig_list.h"
#include "flags.h"
#include "form_ode.h"
#include "ggets.h"
#include "graphics.h"
#include "integrate.h"
#include "load_eqn.h"
#include "main.h"
#include "markov.h"
#include "menudrive.h"
#include "my_rhs.h"
#include "numerics.h"
#include "util/matrixalgebra.h"
#include "util/timeutil.h"

/* --- Forward declarations --- */
static void pr_evec(double *x, double *ev, int n, int pr, double eval,int type);


/* --- Data --- */
static int gear_pivot[MAXODE];
static int ShootType[8];
static double pertst[7][2][3]={{{2,3,1},{2,12,1}},
						{{4.5,6,1},{12,24,1}},
						{{7.333,9.167,.5},{24,37.89,2}},
						{{10.42,12.5,.1667},{37.89,53.33,1}},
						{{13.7,15.98,.04133},{53.33,70.08,.3157}},
						{{17.15,1,.008267},{70.08,87.97,.07407}},
						{{1,1,1},{87.97,1,.0139}}};

int ShootICFlag;
int ShootIndex;
int StableManifoldColor=8;
int UnstableManifoldColor=5;

double ShootIC[8][MAXODE];


/* --- Functions --- */

/* main fixed point finder */
void do_sing(double *x, double eps, double err, double big, int maxit, int n, int *ierr, float *stabinfo) {
	int kmem,i,j,ipivot[MAXODE];
	int oldcol,dummy;
	int rp=0,rn=0,cp=0,cn=0,im=0;
	int pose=0,nege=0,pr;
	double *work,*eval,*b,*bp,*oldwork,*ework;
	double temp,oldt=DELTA_T,old_x[MAXODE];

	char ch;
	double real,imag;
	double bigpos=-1e10,bigneg=1e10;
	int bpos=0,bneg=0;
	/* float xl[MAXODE]; */
	kmem=n*(2*n+5)+50;
	work=(double *)malloc(sizeof(double)*kmem);
	if(work==NULL) {
		err_msg("Insufficient core ");
		return;
	}
	ShootICFlag=0;
	ShootIndex=0;
	for(i=0;i<n;i++) {
		old_x[i]=x[i];
	}
	oldwork=work+n*n;
	eval=oldwork+n*n;
	b=eval+2*n;
	bp=b+n;
	ework=bp+n;
	rooter(x,err,eps,big,work,ierr,maxit,n);
	if(*ierr!=0) {
		free(work);
		err_msg("Could not converge to root");
		for(i=0;i<n;i++) {
			x[i]=old_x[i];
		}
		return;
	}
	ping();
	for(i=0;i<n*n;i++) {
		oldwork[i]=work[i];
	}
	/* Transpose for Eigen        */
	for(i=0;i<n;i++) {
		for(j=i+1;j<n;j++) {
			temp=work[i+j*n];
			work[i+j*n]=work[i*n+j];
			work[i*n+j]=temp;
		}
	}
	eigen(n,work,eval,ework,ierr);
	if(*ierr!=0) {
		err_msg("Could not compute eigenvalues");
		free(work);
		return;
	}
	/* succesfully computed evals now lets work with them */
	ch='n';
	if(!PAR_FOL) {
		ch=(char)TwoChoice("YES","NO","Print eigenvalues?","yn");
	}
	pr=0;

	if(ch=='y') {
		plintf("\n Eigenvalues:\n");
		pr=1;
	}
	for(i=0;i<n;i++) {
		real=eval[2*i];
		imag=eval[2*i+1];
		if(pr==1) {
			plintf(" %f  +  i  %f \n",real,imag);
		}
		if(METHOD==0) {
			real=real*real+imag*imag-1.00;
		}
		if(fabs(imag)<.00000001) {
			imag=0.0;
		}
		if(real<0.0) {
			if(imag!=0.0) {
				cn++;
				if(real<bigneg) {
					bigneg=real;
					bneg=-1;
				}
			} else {
				rn++;
				nege=i;
				if(real<bigneg) {
					bigneg=real;
					bneg=i;
				}
			}
		}
		if(real>0.0) {
			if(imag!=0.0) {
				cp++;
				if(real>bigpos) {
					bigpos=real;
					bpos=-1;
				}
			} else {
				rp++;
				pose=i;
				if(real>bigpos) {
					bigpos=real;
					bpos=i;
				}
			}
		}
		if((real==0.0) && (imag!=0.0)) {
			im++;
		}
	}     /* eigenvalue count */
	if(((rp+cp)!=0) && ((rn+cn)!=0)) {
		eq_symb(x,1);
	} else {
		if((rp+cp)!=0) {
			eq_symb(x,0);
		} else {
			eq_symb(x,3);
		}
	}
	*stabinfo=(float)(cp+rp)+(float)(cn+rn)/1000.0;

	/* Lets change Work back to transposed oldwork */
	for(i=0;i<n;i++) {
		for(j=i+1;j<n;j++) {
			temp=oldwork[i+j*n];
			work[i+j*n]=oldwork[i*n+j];
			work[i*n+j]=temp;
		}
	}
	create_eq_box(cp,cn,rp,rn,im,x,eval,n);
	if(((rp==1) || (rn==1)) && (n>1)) {
		ch='n';
		if(!PAR_FOL) {
			ch=(char)TwoChoice("YES","NO","Draw Invariant Sets?","yn");
		}
		if((ch=='y') || (PAR_FOL&&SHOOT)) {
			oldt=DELTA_T;
			if(rp==1) {
				/* plintf(" One real positive -- pos=%d lam=%g \n",pose,eval[2*pose]); */
				/*     for(i=0;i<n*n;i++)printf(" w=%g o=%g \n",work[i],oldwork[i]); */
				get_evec(work,oldwork,b,bp,n,maxit,err,ipivot,eval[2*pose],ierr);
				if(*ierr==0) {
					change_current_linestyle(UnstableManifoldColor,&oldcol);
					pr_evec(x,b,n,pr,eval[2*pose],1);
					DELTA_T=fabs(DELTA_T);
					shoot(bp,x,b,1);
					shoot(bp,x,b,-1);
					change_current_linestyle(oldcol,&dummy);
				} else {
					err_msg("Failed to compute eigenvector");
				}
			}
			if(rn==1) {
				get_evec(work,oldwork,b,bp,n,maxit,err,ipivot,eval[2*nege],ierr);
				if(*ierr==0) {
					change_current_linestyle(StableManifoldColor,&oldcol);
					pr_evec(x,b,n,pr,eval[2*nege],-1);
					DELTA_T=-fabs(DELTA_T);
					shoot(bp,x,b,1);
					shoot(bp,x,b,-1);
					change_current_linestyle(oldcol,&dummy);
				} else {
					err_msg("Failed to compute eigenvector");
				}
			}
			DELTA_T=oldt;
		}
	}  /* end of normal shooting stuff */

	/* strong (un) stable manifold calculation
	only one-d manifolds calculated */
	/* lets check to see if it is relevant */
	if(((rn>1) && (bneg>=0)) || ((rp>1) && (bpos>=0))) {
		ch='n';
		if(!PAR_FOL) {
			ch=(char)TwoChoice("YES","NO","Draw Strong Sets?","yn");
		}
		if((ch=='y')||(PAR_FOL&&SHOOT)) {
			oldt=DELTA_T;
			if((rp>1) && (bpos>=0)) {/* then there is a strong unstable */
				plintf("strong unstable %g \n",bigpos);
				get_evec(work,oldwork,b,bp,n,maxit,err,ipivot,bigpos,ierr);
				if(*ierr==0) {
					change_current_linestyle(UnstableManifoldColor,&oldcol);
					pr_evec(x,b,n,pr,bigpos,1);
					DELTA_T=fabs(DELTA_T);
					shoot(bp,x,b,1);
					shoot(bp,x,b,-1);
					change_current_linestyle(oldcol,&dummy);
				} else {
					err_msg("Failed to compute eigenvector");
				}
			}
			if((rn>1) && (bneg>=0)) {/* then there is a strong stable */
				plintf("strong stable %g \n",bigneg);
				get_evec(work,oldwork,b,bp,n,maxit,err,ipivot,bigneg,ierr);
				if(*ierr==0) {
					change_current_linestyle(StableManifoldColor,&oldcol);
					pr_evec(x,b,n,pr,bigneg,-1);
					DELTA_T=-fabs(DELTA_T);
					shoot(bp,x,b,1);
					shoot(bp,x,b,-1);
					change_current_linestyle(oldcol,&dummy);
				} else {
					err_msg("Failed to compute eigenvector");
				}
			}
		}
		DELTA_T=oldt;
	}
	free(work);
	return;
}


/* fixed point with no requests and store manifolds */
void do_sing_info(double *x, double eps, double err, double big,
				  int maxit, int n, double *er, double *em, int *ierr) {
	int kmem,i,j,ipivot[MAXODE];

	int rp=0,rn=0,cp=0,cn=0,im=0;
	int pose=0,nege=0,pr=0;
	double *work,*eval,*b,*bp,*oldwork,*ework;
	double temp,old_x[MAXODE];
	double real,imag;
	double bigpos=-1e10,bigneg=1e10;

	/* float xl[MAXODE]; */
	kmem=n*(2*n+5)+50;
	work=(double *)malloc(sizeof(double)*kmem);
	if(work==NULL) {
		/* err_msg("Insufficient core "); */
		return;
	}
	ShootICFlag=0;
	ShootIndex=0;
	for(i=0;i<n;i++) {
		old_x[i]=x[i];
	}
	oldwork=work+n*n;
	eval=oldwork+n*n;
	b=eval+2*n;
	bp=b+n;
	ework=bp+n;
	rooter(x,err,eps,big,work,ierr,maxit,n);
	if(*ierr!=0) {
		free(work);
		/* err_msg("Could not converge to root"); */
		for(i=0;i<n;i++) {
			x[i]=old_x[i];
		}
		return;
	}

	for(i=0;i<n*n;i++) {
		oldwork[i]=work[i];
	}
	/* Transpose for Eigen        */
	for(i=0;i<n;i++) {
		for(j=i+1;j<n;j++) {
			temp=work[i+j*n];
			work[i+j*n]=work[i*n+j];
			work[i*n+j]=temp;
		}
	}
	eigen(n,work,eval,ework,ierr);
	if(*ierr!=0) {
		free(work);
		return;
	}
	/* succesfully computed evals now lets work with them */

	for(i=0;i<n;i++) {
		real=eval[2*i];
		imag=eval[2*i+1];
		er[i]=real;
		em[i]=imag;

		if(METHOD==0) {
			real=real*real+imag*imag-1.00;
		}
		if(fabs(imag)<.00000001) {
			imag=0.0;
		}
		if(real<0.0) {
			if(imag!=0.0) {
				cn++;
				if(real<bigneg) {
					bigneg=real;
				}
			} else {
				rn++;
				nege=i;
				if(real<bigneg) {
					bigneg=real;
				}
			}
		}
		if(real>0.0) {
			if(imag!=0.0) {
				cp++;
				if(real>bigpos) {
					bigpos=real;
				}
			} else {
				rp++;
				pose=i;
				if(real>bigpos) {
					bigpos=real;
				}
			}
		}
		if((real==0.0) && (imag!=0.0)) {
			im++;
		}
	}     /* eigenvalue count */
	if(((rp+cp)!=0) && ((rn+cn)!=0)) {
		eq_symb(x,1);
	} else {
		if((rp+cp)!=0) {
			eq_symb(x,0);
		} else {eq_symb(x,3);
		}
	}
	/* Lets change Work back to transposed oldwork */
	for(i=0;i<n;i++) {
		for(j=i+1;j<n;j++) {
			temp=oldwork[i+j*n];
			work[i+j*n]=oldwork[i*n+j];
			work[i*n+j]=temp;
		}
	}
	if((n>1)) {
		if(rp==1) {
			/* plintf(" One real positive -- pos=%d lam=%g \n",pose,eval[2*pose]); */
			/*     for(i=0;i<n*n;i++)printf(" w=%g o=%g \n",work[i],oldwork[i]); */
			get_evec(work,oldwork,b,bp,n,maxit,err,ipivot,eval[2*pose],ierr);

			if(*ierr==0) {
				pr_evec(x,b,n,pr,eval[2*pose],1);
			}
		}

		if(rn==1) {
			get_evec(work,oldwork,b,bp,n,maxit,err,ipivot,eval[2*nege],ierr);
			if(*ierr==0) {
				pr_evec(x,b,n,pr,eval[2*nege],-1);
			}
		}
	}

	free(work);
	return;
}


int gear(int n, double *t, double tout, double *y, double hmin, double hmax,
		 double eps, int mf, double *error, int *kflag, int *jstart,
		 double *work, int *iwork) {
	if(NFlags==0) {
		return(ggear( n,t, tout,y, hmin, hmax,eps,
					  mf,error,kflag,jstart,work,iwork));
	}
	return(one_flag_step_gear(n,t, tout,y, hmin,
							  hmax,eps,mf,error,kflag,jstart,work,iwork));
}


int ggear(int n, double *t, double tout, double *y, double hmin, double hmax,
		  double eps, int mf, double *error, int *kflag, int *jstart,
		  double *work, int *iwork) {
	double deltat=0.0,hnew=0.0,hold=0.0,h=0.0,racum=0.0,told=0.0,r=0.0,d=0.0;
	double *a,pr1,pr2,pr3,r1;
	double *dermat,*save[8],*save9,*save10,*save11,*save12;
	double enq1=0.0,enq2=0.0,enq3=0.0,pepsh=0.0,e=0.0,edwn=0.0,eup=0.0,bnd=0.0;
	double *ytable[8],*ymax,*work2;
	int i,iret=0,maxder=0,j=0,k=0,iret1=0,nqold=0,nq=0,newq=0;
	int idoub=0,mtyp=0,iweval=0,j1=0,j2=0,l=0,info=0,job=0,nt=0;

	for(i=0;i<8;i++) {
		save[i]=work+i*n;
		ytable[i]=work+(8+i)*n;
	}
	save9=work+16*n;
	save10=work+17*n;
	save11=work+18*n;
	save12=work+19*n;
	ymax=work+20*n;
	dermat=work+21*n;
	a=work+21*n+n*n;
	work2=work+21*n+n*n+10;
	if(*jstart!=0) {
		k=iwork[0];
		nq=iwork[1];
		nqold=iwork[2];
		idoub=iwork[3];
		maxder=6;
		mtyp=1;
		iret1=iwork[4];
		iret=iwork[5];
		newq=iwork[6];
		iweval=iwork[7];
		hold=work2[1];
		h=work2[0];
		hnew=work2[2];
		told=work2[3];
		racum=work2[4];
		enq1=work2[5];
		enq2=work2[6];
		enq3=work2[7];
		pepsh=work2[8];
		e=work2[9];
		edwn=work2[10];
		eup=work2[11];
		bnd=work2[12];
	}
	deltat=tout-*t;
	if(*jstart==0) {
		h=SIGN(hmin,deltat);
	}
	if(fabs(deltat)<hmin) {
		return(-1);
	}
	maxder=6;
	for(i=0;i<n;i++) {
		ytable[0][i]=y[i];
	}
L70:
	iret=1;
	*kflag=1;
	if((h>0.0) && ((*t+h)>tout)) {
		h=tout-*t;
	}
	if((h<0.0) && ((*t+h)<tout)) {
		h=tout-*t;
	}
	if(*jstart<=0) {
		goto L120;
	}
L80:
	for(i=0;i<n;i++) {
		for(j=1;j<=k;j++) {
			save[j-1][i]=ytable[j-1][i];
		}
	}
	hold=hnew;
	if(h==hold) {
		goto L110;
	}
L100:
	racum=h/hold;
	iret1=1;
	goto L820;
L110:
	nqold=nq;
	told= *t;
	racum=1.0;
	if(*jstart>0) {
		goto L330;
	}
	goto L150;
L120:
	if(*jstart==-1) {
		goto L140;
	}
	nq=1;
	rhs(*t,ytable[0],save11,n);
	for(i=0;i<n;i++) {
		ytable[1][i]=save11[i]*h;
		ymax[i]=1.00;
	}
	hnew=h;
	k=2;
	goto L80;
L140:
	if(nq==nqold) {
		*jstart=1;
	}
	*t=told;
	nq=nqold;
	k=nq+1;
	goto L100;
L150:
	if(nq>6) {
		*kflag=-2;
		goto L860;
	}
	switch(nq) {
	case 1:
		a[0]=-1.00;
		a[1]=-1.00;
		break;
	case 2:
		a[0]=-2.0/3.0;
		a[1]=-1.00;
		a[2]=-1.0/3.0;
		break;
	case 3:
		a[0]=-6.0/11.0;
		a[1]=-1.00;
		a[2]=a[0];
		a[3]=-1.0/11.0;
		break;
	case 4:
		a[0]=-.48;
		a[1]=-1.00;
		a[2]=-.70;
		a[3]=-.2;
		a[4]=-.02;
		break;
	case 5:
		a[0]=-120.0/274.;
		a[1]=-1.00;
		a[2]=-225./274.;
		a[3]=-85./274.;
		a[4]=-15./274.;
		a[5]=-1./274.;
		break;
	case 6:
		a[0]=-180./441.;
		a[1]=-1.0;
		a[2]=-58./63.;
		a[3]=-5./12.;
		a[4]=-25./252.;
		a[5]=-3./252.;
		a[6]=-1./1764;
		break;
	}
	/*L310:*/
	k=nq+1;
	idoub=k;
	mtyp=(4-mf)/2;
	enq2=.5/(double)(nq+1);
	enq3=.5/(double)(nq+2);
	enq1=.5/(double)nq;
	pepsh=eps;
	eup=sqr2(pertst[nq-1][0][1]*pepsh);
	e=sqr2(pertst[nq-1][0][0]*pepsh);
	edwn=sqr2(pertst[nq-1][0][2]*pepsh);
	if(edwn==0.0) {
		goto L850;
	}
	bnd=eps*enq3/(double)n;

	/*L320:*/
	iweval=2;
	if(iret==2) {
		goto L750;
	}
L330:
	*t=*t+h;
	for(j=2;j<=k;j++) {
		for(j1=j;j1<=k;j1++) {
			j2=k-j1+j-1;
			for(i=0;i<n;i++) {
				ytable[j2-1][i]=ytable[j2-1][i]+ytable[j2][i];
			}
		}
	}
	for(i=0;i<n;i++) {
		error[i]=0.0;
	}
	for(l=0;l<3;l++) {
		rhs(*t,ytable[0],save11,n);
		if(iweval<1) {
			/*  plintf("iweval=%d \n",iweval);
		 for(i=0;i<n;i++)printf("up piv = %d \n",gear_pivot[i]);*/
			goto L460;
		}
		/*       JACOBIAN COMPUTED   */
		for(i=0;i<n;i++) {
			save9[i]=ytable[0][i];
		}
		for(j=0;j<n;j++) {
			r=eps*MAX(eps,fabs(save9[j]));
			ytable[0][j]=ytable[0][j]+r;
			d=a[0]*h/r;
			rhs(*t,ytable[0],save12,n);
			for(i=0;i<n;i++) {
				dermat[n*i+j]=(save12[i]-save11[i])*d;
			}
			ytable[0][j]=save9[j];

		}
		for(i=0;i<n;i++) {
			dermat[n*i+i]+=1.0;
		}
		iweval=-1;
		sgefa(dermat,n,n,gear_pivot,&info);

		if(info==-1) {
			j1=1;
		} else {
			j1=-1;
		}
		if(j1<0) {
			goto L520;
		}
L460:
		for(i=0;i<n;i++) {
			save12[i]=ytable[1][i]-save11[i]*h;
		}
		for(i=0;i<n;i++) {
			save9[i]=save12[i];
		}
		job=0;
		sgesl(dermat,n,n,gear_pivot,save9,job);
		nt=n;
		for(i=0;i<n;i++) {
			ytable[0][i]=ytable[0][i]+a[0]*save9[i];
			ytable[1][i]=ytable[1][i]-save9[i];
			error[i]+=save9[i];
			if(fabs(save9[i])<=(bnd*ymax[i])) {
				nt--;
			}
		}
		if(nt<=0) {
			goto L560;
		}
	}

L520:
	/*        UH Oh */
	*t=told;
	if((h<=(hmin*1.000001)) && ((iweval-mtyp)<-1)) {
		goto L530;
	}
	if(iweval!=0) {
		racum*=.25;
	}
	iweval=mf;
	iret1=2;
	goto L820;

L530:
	*kflag=-3;

L540:
	for(i=0;i<n;i++) {
		for(j=1;j<=k;j++) {
			ytable[j-1][i]=save[j-1][i];
		}
	}
	h=hold;
	nq=nqold;
	*jstart=nq;
	goto L860;

L560:
	d=0.0;
	for(i=0;i<n;i++) {
		d+=sqr2(error[i]/ymax[i]);
	}
	iweval=0;
	if(d>e) {
		goto L610;
	}
	if(k>=3) {
		for(j=3;j<=k;j++) {
			for(i=0;i<n;i++) {
				ytable[j-1][i]=ytable[j-1][i]+a[j-1]*error[i];
			}
		}
	}
	*kflag=1;
	hnew=h;
	if(idoub<=1) {
		goto L620;
	}
	idoub--;
	if(idoub<=1) {
		for(i=0;i<n;i++) {
			save10[i]=error[i];
		}
	}
	goto L770;

L610:
	*kflag-=2;
	if(h<=hmin*1.00001) {
		goto L810;
	}
	*t=told;
	if(*kflag<=-5) {
		goto L790;
	}

L620:
	pr2=1.2*pow(d/e,enq2);
	pr3=1.0e20;
	if((nq<maxder)&&(*kflag>-1)) {
		d=0.0;
		for(i=0;i<n;i++) {
			d+=sqr2((error[i]-save10[i])/ymax[i]);
		}
		pr3=1.4*pow(d/eup,enq3);
	}
	pr1=1.0e20;
	if(nq>1) {
		d=0.0;
		for(i=0;i<n;i++) {
			d+=sqr2(ytable[k-1][i]/ymax[i]);
		}
		pr1=1.3*pow(d/edwn,enq1);
	}
	if(pr2<=pr3) {
		goto L720;
	}
	if(pr3<pr1) {
		goto L730;
	}
L670:
	r=1.0/MAX(pr1,0.0001);
	newq=nq-1;

L680:
	idoub=10;
	if((*kflag==1)&&(r<1.1)) {
		goto L770;
	}
	if(newq<=nq) {
		goto L700;
	}
	for(i=0;i<n;i++) {
		ytable[newq][i]=error[i]*a[k-1]/(double)k;
	}

L700:
	k=newq+1;
	if(*kflag==1) {
		goto L740;
	}
	racum=racum*r;
	iret1=3;
	goto L820;

L710:
	if(newq==nq) {
		goto L330;
	}
	nq=newq;
	goto L150;

L720:
	if(pr2>pr1) {
		goto L670;
	}
	newq=nq;
	r=1.0/MAX(pr2,.0001);
	goto L680;

L730:
	r=1.0/MAX(pr3,.0001);
	newq=nq+1;
	goto L680;

L740:
	iret=2;
	h=h*r;
	hnew=h;
	if(nq==newq) {
		goto L750;
	}
	nq=newq;
	goto L150;

L750:
	r1=1.0;
	for(j=2;j<=k;j++) {
		r1=r1*r;
		for(i=0;i<n;i++) {
			ytable[j-1][i]=ytable[j-1][i]*r1;
		}
	}
	idoub=k;

L770:
	for(i=0;i<n;i++) {
		ymax[i]=MAX(ymax[i],fabs(ytable[0][i]));
	}
	*jstart=nq;
	if((h>0.0)&&(*t>=tout)) {
		goto L860;
	}
	if((h<0.0)&&(*t<=tout)) {
		goto L860;
	}
	goto L70;

L790:

	if(nq==1) {
		goto L850;
	}
	rhs(*t,ytable[0],save11,n);
	r=h/hold;
	for(i=0;i<n;i++) {
		ytable[0][i]=save[0][i];
		save[1][i]=hold*save11[i];
		ytable[1][i]=r*save[1][i];
	}
	nq=1;
	*kflag=1;
	goto L150;

L810:
	*kflag=-1;
	hnew=h;
	*jstart=nq;
	goto L860;

L820:
	racum=MAX(fabs(hmin/hold),racum);
	racum=MIN(racum,fabs(hmax/hold));
	r1=1.0;
	for(j=2;j<=k;j++) {
		r1=r1*racum;
		for(i=0;i<n;i++) {
			ytable[j-1][i]=save[j-1][i]*r1;
		}
	}
	h=hold*racum;
	for(i=0;i<n;i++) {
		ytable[0][i]=save[0][i];
	}
	idoub=k;
	if(iret1==1) {
		goto L110;
	}
	if(iret1==2) {
		goto L330;
	}
	if(iret1==3) {
		goto L710;
	}

L850:
	*kflag=-4;
	goto L540;

L860:
	for(i=0;i<n;i++) {
		y[i]=ytable[0][i];
	}
	iwork[0]=k;
	iwork[1]=nq;
	iwork[2]=nqold;
	work2[0]=h;
	work2[1]=hold;
	work2[2]=hnew;
	work2[3]=told;
	work2[4]=racum;
	work2[5]=enq1;
	work2[6]=enq2;
	work2[7]=enq3;
	work2[8]=pepsh;
	work2[9]=e;
	work2[10]=edwn;
	work2[11]=eup;
	work2[12]=bnd;
	iwork[3]=idoub;
	iwork[4]=iret1;
	iwork[5]=iret;
	iwork[6]=newq;
	iwork[7]=iweval;
	return(1);
}


void rooter(double *x, double err, double eps, double big, double *work, int *ierr, int maxit, int n) {
	int i,iter,ipivot[MAXODE],info;
	char ch;
	double *xp,*yp,*y,*xg,*dermat,*dely;
	double r;
	dermat=work;
	xg=dermat+n*n;
	yp=xg+n;
	xp=yp+n;
	y=xp+n;
	dely=y+n;
	iter=0;
	*ierr=0;
	while(1) {
		ch=my_abort();

		if(ch==27) {
			*ierr=1;
			return;
		}
		if(ch=='/') {
			*ierr=1;
			ENDSING=1;
			return;
		}
		if(ch=='p') {
			PAUSER=1;
		}
		getjac(x,y,yp,xp,eps,dermat,n);
		sgefa(dermat,n,n,ipivot,&info);
		if(info!=-1) {
			*ierr=1;
			return;
		}
		for(i=0;i<n;i++) {
			dely[i]=y[i];
		}
		sgesl(dermat,n,n,ipivot,dely,0);
		r=0.0;
		for(i=0;i<n;i++) {
			x[i]=x[i]-dely[i];
			r=r+fabs(dely[i]);
		}
		if(r<err) {
			getjac(x,y,yp,xp,eps,dermat,n);
			if(METHOD==0) {
				for(i=0;i<n;i++) {
					dermat[i*(n+1)]+=1.0;
				}
			}
			return; /* success !! */
		}
		if((r/(double)n)>big) {
			*ierr=1;
			return;
		}
		iter++;
		if(iter>maxit) {
			*ierr=1;
			return;
		}
	}
}


/* this uses the current labeled saddle point stuff to integrate */
void shoot_this_now(void) {
	int i,k,type,oldcol,dummy;
	double x[MAXODE],olddt;
	if(ShootIndex<1)
		return;

	olddt=DELTA_T;

	for(k=0;k<ShootIndex;k++) {
		for(i=0;i<NODE;i++) {
			x[i]=ShootIC[k][i];
		}
		type=ShootType[k];
		if(type>0) {
			change_current_linestyle(UnstableManifoldColor,&oldcol);
			DELTA_T=fabs(DELTA_T);
			shoot_easy(x);
			change_current_linestyle(oldcol,&dummy);
		}
		if(type<0) {
			change_current_linestyle(StableManifoldColor,&oldcol);
			DELTA_T=-fabs(DELTA_T);
			shoot_easy(x);
			change_current_linestyle(oldcol,&dummy);
		}
	}
	DELTA_T=olddt;
}


/* --- Static functions --- */
static void pr_evec(double *x, double *ev, int n, int pr, double eval,int type) {
	int i;
	double d=fabs(DELTA_T)*.1;
	ShootICFlag=1;
	if(ShootIndex<7) {
		for(i=0;i<n;i++) {
			ShootIC[ShootIndex][i]=x[i]+d*ev[i];
			ShootType[ShootIndex]=type;
			ShootIC[ShootIndex+1][i]=x[i]-d*ev[i];
			ShootType[ShootIndex+1]=type;
		}
		ShootIndex+=2;
	}
	if(pr==0) {
		return;
	}
}
