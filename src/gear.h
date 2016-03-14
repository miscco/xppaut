#ifndef XPPAUT_GEAR_H
#define XPPAUT_GEAR_H

#include "xpplim.h"

/* --- Data --- */
extern double ShootIC[8][MAXODE];

extern int ShootICFlag;
extern int ShootIndex;
extern int StableManifoldColor;
extern int UnstableManifoldColor;

/* --- Functions --- */
void do_sing(double *x, double eps, double err, double big, int maxit, int n, int *ierr, float *stabinfo);
void do_sing_info(double *x, double eps, double err, double big, int maxit, int n, double *er, double *em, int *ierr);
int gear(int n, double *t, double tout, double *y, double hmin, double hmax, double eps, int mf, double *error, int *kflag, int *jstart, double *work, int *iwork);
int ggear(int n, double *t, double tout, double *y, double hmin, double hmax, double eps, int mf, double *error, int *kflag, int *jstart, double *work, int *iwork);
void rooter(double *x, double err, double eps, double big, double *work, int *ierr, int maxit, int n);
void shoot_this_now(void);

#endif /* XPPAUT_GEAR_H */
