#ifndef XPPAUT_FLAGS_H
#define XPPAUT_FLAGS_H

/* --- Data --- */
extern int NFlags;
extern double STOL;

/* --- Functions --- */
int add_global(char *cond, int sign, char *rest);
void show_flags(void);
int compile_flags(void);
int one_flag_step(double *yold, double *ynew, int *istart, double told, double *tnew, int neq, double *s);
int one_flag_step_adap(double *y, int neq, double *t, double tout, double eps, double *hguess, double hmin, double *work, int *ier, double epjac, int iflag, int *jstart);

#endif /* XPPAUT_FLAGS_H */

