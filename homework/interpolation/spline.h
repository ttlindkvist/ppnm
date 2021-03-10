#ifndef __SPLINE_H__
#define __SPLINE_H__

double linterp(int n, double *x, double *y, double z);
double linterp_integ(int n, double *x, double *y, double z);

typedef struct {int n; double *x, *y, *b, *c;} qspline;
qspline *qspline_init(int n, double *x, double *y);
void qspline_free(qspline *s);
double qspline_eval(qspline *s, double z);
double qspline_deriv(qspline *s, double z);
double qspline_integ(qspline *s, double z);

#endif