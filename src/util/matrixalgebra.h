#ifndef XPPAUT_UTIL_MATRIXALGEBRA_H
#define XPPAUT_UTIL_MATRIXALGEBRA_H

/* --- Functions --- */
int bandfac(double *a, int ml, int mr, int n);
void bandsol(double *a, double *b, int ml, int mr, int n);
void get_the_jac(double t, double *y, double *yp, double *ypnew, double *dfdy, int neq, double eps, double scal);

#endif // XPPAUT_UTIL_MATRIXALGEBRA_H
