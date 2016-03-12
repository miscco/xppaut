#include "matrixalgebra.h"

#include <math.h>

#include "../main.h"
#include "../numerics.h"

/* --- Macros --- */
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))

/* --- Forward declarations --- */
static void get_band_jac(double *a, double *y, double t, double *ypnew, double *ypold, int n, double eps, double scal);


/* --- Functions --- */
/* factors the matrix */
int bandfac(double *a, int ml, int mr, int n) {
	int i,j,k;
	int row,rowi,m,r0,ri0;
	double al;

	int n1=n-1;
	int mt=ml+mr+1;

	for(row=0;row<n;row++) {
		r0=row*mt+ml;
		if((al=a[r0])==0.0) {
			return(-1-row);
		}
		al=1.0/al;
		m=MIN(mr,n1-row);
		for(j=1;j<=m;j++) {
			a[r0+j]=a[r0+j]*al;
		}
		a[r0]=al;
		for(i=1;i<=ml;i++) {
			rowi=row+i;
			if(rowi>n1) {
				break;
			}
			ri0=rowi*mt+ml;
			al=a[ri0-i];
			if(al==0.0) {
				continue;
			}
			for(k=1;k<=m;k++) {
				a[ri0-i+k]=a[ri0-i+k]-(al*a[r0+k]);
			}
			a[ri0-i]=-al;
		}
	}
	return(0);
}


/* requires that the matrix be factored */
void bandsol(double *a, double *b, int ml, int mr, int n) {
	int i,j,k,r0;

	int mt=ml+mr+1;
	int m,n1=n-1,row;

	for(i=0;i<n;i++) {
		r0=i*mt+ml;
		m=MAX(-ml,-i);
		for(j=m;j<0;j++) {
			b[i] += a[r0+j]*b[i+j];
		}
		b[i] *= a[r0];
	}
	for(row=n1-1;row>=0;row--) {
		m=MIN(mr,n1-row);
		r0=row*mt+ml;
		for(k=1;k<=m;k++) {
			b[row]=b[row]-a[r0+k]*b[row+k];
		}
	}
}


/* this assumes that yp is already computed */
void get_the_jac(double t, double *y, double *yp, double *ypnew, double *dfdy, int neq, double eps, double scal) {
	int i,j;
	double yold,del,dsy;
	if(cv_bandflag) {
		get_band_jac(dfdy,y,t,ypnew,yp,neq,eps,scal);
	} else {
		for(i=0;i<neq;i++) {
			del=eps*MAX(eps,fabs(y[i]));
			dsy=scal/del;
			yold=y[i];
			y[i]=y[i]+del;
			rhs(t,y,ypnew,neq);
			for(j=0;j<neq;j++) {
				dfdy[j*neq+i]=dsy*(ypnew[j]-yp[j]);
			}
			y[i]=yold;
		}
	}
}


/* --- Static functions --- */
static void get_band_jac(double *a, double *y, double t, double *ypnew, double *ypold, int n, double eps, double scal) {
	int ml=cv_bandlower,mr=cv_bandupper;
	int i,j,k,n1=n-1,mt=ml+mr+1;
	double yhat;
	double dy;
	double dsy;

	for(i=0;i<(n*mt);i++) {
		a[i]=0.0;
	}
	for(i=0;i<n;i++) {
		yhat=y[i];
		dy=eps*(eps+fabs(yhat));
		dsy=scal/dy;
		y[i] += dy;
		rhs(t,y,ypnew,n);
		for(j=-ml;j<=mr;j++) {
			k=i-j;
			if(k<0 || k>n1) {
				continue;
			}
			a[k*mt+j+ml]=dsy*(ypnew[k]-ypold[k]);
		}
		y[i]=yhat;
	}
}
