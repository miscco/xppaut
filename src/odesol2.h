#ifndef XPPAUT_ODESOL2_H
#define XPPAUT_ODESOL2_H

/* --- Functions --- */
int adams(double *y, double *tim, double dt, int nstep, int neq, int *ist, double *work);
int rb23(double *y, double *tstart, double tfinal, int *istart, int n, double *work, int *ierr);
int rosen(double *y, double *tstart, double tfinal, int *istart, int n, double *work, int *ierr);

#endif /* XPPAUT_ODESOL2_H */
