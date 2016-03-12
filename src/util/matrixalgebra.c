#include "matrixalgebra.h"

#include <math.h>

#include "../main.h"
#include "../numerics.h"

/* --- Macros --- */
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))

/* --- Forward declarations --- */
static void get_band_jac(double *a, double *y, double t, double *ypnew, double *ypold, int n, double eps, double scal);
static int isamax(int n, double *sx, int incx);
static void saxpy(int n, double sa, double *sx, int incx, double *sy, int incy);
static double sdot(int n, double *sx, int incx, double *sy, int incy);
static void sscal(int n, double sa, double *sx, int incx);


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


/**
 * Factors a real matrix by Gaussian elimination.
 *
 * @param a (in/out) the matrix to be factored.
 * @param lda the leading dimension of the array A.
 * @param n the order of the matrix A.
 * @param ipvt (out) an integer vector of pivot indices.
 * @param info (out) = 0  normal value.
 *             = k  if  u(k,k) .eq. 0.0 .  This is not an error
 *                  condition for this subroutine, but it does
 *                  indicate that SGESL or SGEDI will divide by zero
 *                  if called.  Use  RCOND  in SGECO for a reliable
 *                  indication of singularity.
 */
void sgefa(double *a, int lda, int n, int *ipvt, int *info) {
	int j,k,kp1,l,nm1;
	double t;
	*info=-1;
	nm1=n-1;
	if(nm1>0) {
		for(k=1;k<=nm1;k++) {
			kp1=k+1;
			l=isamax(n-k+1,&a[(k-1)*lda+k-1],lda)+k-1;
			ipvt[k-1]=l;
			if(a[l*lda+k-1]!=0.0) {
				if(l!=(k-1)) {
					t=a[l*lda+k-1];
					a[l*lda+k-1]=a[(k-1)*lda+k-1];
					a[(k-1)*lda+k-1]=t;
				}
				t=-1.0/a[(k-1)*lda+k-1];
				sscal(n-k,t,(a+k*lda+k-1),lda);
				for(j=kp1;j<=n;j++)	{
					t=a[l*lda+j-1];
					if(l!=(k-1)) {
						a[l*lda+j-1]=a[(k-1)*lda+j-1];
						a[(k-1)*lda+j-1]=t;
					}
					saxpy(n-k,t,(a+k*lda+k-1),lda,(a+k*lda+j-1),lda);
				}
			} else {
				*info=k-1;
			}
		}
	}
	ipvt[n-1]=n-1;
	if(a[(n-1)*lda+n-1]==0.0) {
		*info=n-1;
	}
}


/**
 * Solves the real system
 *   A * X = B  or  transpose(A) * X = B
 * using factors computed by sgeco or sgefa.
 *
 * @param a the output from sgeco or sgefa.
 * @param lda the leading dimension of the array A.
 * @param n the order of the matrix A.
 * @param ipvt the pivot vector from sgeco or sgefa.
 * @param b (in/out) the right hand side vector.
 * @param job whether to solve the transpose(A) problem.
 */
void sgesl(double *a, int lda, int n, int *ipvt, double *b, int job) {
	int k,kb,l,nm1;
	double t;
	nm1=n-1;
	if(job==0) {
		if(nm1>=1) {
			for(k=1;k<=nm1;k++) {
				l=ipvt[k-1];
				t=b[l];
				if(l!=(k-1)) {
					b[l]=b[k-1];
					b[k-1]=t;
				}
				saxpy(n-k,t,(a+lda*k+k-1),lda,(b+k),1);
			}
		}
		for(kb=1;kb<=n;kb++) {
			k=n+1-kb;
			b[k-1]=b[k-1]/a[(k-1)*lda+k-1];
			t=-b[k-1];
			saxpy(k-1,t,(a+k-1),lda,b,1);
		}
		return;
	}
	for(k=1;k<=n;k++) {
		t=sdot(k-1,(a+k-1),lda,b,1);
		b[k-1]=(b[k-1]-t)/a[(k-1)*lda+k-1];
	}
	if(nm1>0) {
		for(kb=1;kb<=nm1;kb++) {
			k=n-kb;
			b[k-1]=b[k-1]+sdot(n-k,(a+k*lda+k-1),lda,b+k,1);
			l=ipvt[k-1];
			if(l!=(k-1)) {
				t=b[l];
				b[l]=b[k-1];
				b[k-1]=t;
			}
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


static int isamax(int n, double *sx, int incx) {
	int i,ix,imax;
	double smax;
	if(n<1)
		return(-1);
	if(n==1)
		return(0);
	if(incx!=1) {
		ix=0;
		imax=0;
		smax=fabs(sx[0]);
		ix+=incx;
		for(i=1;i<n;i++,ix+=incx) {
			if(fabs(sx[ix])>smax) {
				imax=i;
				smax=fabs(sx[ix]);
			}
		}
		return(imax);
	}
	imax=0;
	smax=fabs(sx[0]);
	for(i=1;i<n;i++) {
		if(fabs(sx[i])>smax) {
			imax=i;
			smax=fabs(sx[i]);
		}
	}
	return(imax);
}


static void saxpy(int n, double sa, double *sx, int incx, double *sy, int incy) {
	int i,ix,iy;
	if(n<=0)
		return;
	if(sa==0.0)
		return;
	ix=0;
	iy=0;
	if(incx<0) {
		ix=-n*incx;
	}
	if(incy<0) {
		iy=-n*incy;
	}
	for(i=0;i<n;i++,ix+=incx,iy+=incy) {
		sy[iy]=sy[iy]+sa*sx[ix];
	}
}


static double sdot(int n, double *sx, int incx, double *sy, int incy) {
	int i,ix,iy;
	double stemp=0.0;
	if(n<=0)
		return(0.0);
	ix=0;
	iy=0;
	if(incx<0) {
		ix=-n*incx;
	}
	if(incy<0) {
		iy=-n*incy;
	}
	for(i=0;i<n;i++,ix+=incx,iy+=incy) {
		stemp+=sx[ix]*sy[iy];
	}
	return(stemp);
}


static void sscal(int n, double sa, double *sx, int incx) {
	int i,nincx;
	if(n<=0)
		return;
	nincx=n*incx;
	for(i=0;i<nincx;i+=incx) {
		sx[i]*=sa;
	}
}
