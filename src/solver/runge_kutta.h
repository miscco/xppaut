#ifndef XPPAUT_SOLVER_RUNGE_KUTTA_H
#define XPPAUT_SOLVER_RUNGE_KUTTA_H

int runge_kutta(double *y, double *tim, double dt, int nt, int neq, int *istart, double *work);

#endif // XPPAUT_SOLVER_RUNGE_KUTTA_H
