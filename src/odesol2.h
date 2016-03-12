#ifndef XPPAUT_ODESOL2_H
#define XPPAUT_ODESOL2_H

/* --- Functions --- */
void one_step_euler(double *y, double dt, double *yp, int neq, double *t);
void one_step_rk4(double *y, double dt, double *yval[3], int neq, double *tim);
void one_step_heun(double *y, double dt, double *yval[2], int neq, double *tim);
int euler(double *y, double *tim, double dt, int nt, int neq, int *istart, double *work);
int mod_euler(double *y, double *tim, double dt, int nt, int neq, int *istart, double *work);
int rung_kut(double *y, double *tim, double dt, int nt, int neq, int *istart, double *work);
int adams(double *y, double *tim, double dt, int nstep, int neq, int *ist, double *work);
int rb23(double *y, double *tstart, double tfinal, int *istart, int n, double *work, int *ierr);
int rosen(double *y, double *tstart, double tfinal, int *istart, int n, double *work, int *ierr);

#endif /* XPPAUT_ODESOL2_H */
