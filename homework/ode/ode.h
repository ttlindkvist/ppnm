#ifndef __ODE_H__
#define __ODE_H__
#include<gsl/gsl_matrix.h>

void rkstep45(
	void (*f)(double, gsl_vector*, gsl_vector*), /* the f from dy/dt=f(t,y) */
	double t,              /* the current value of the variable */
	gsl_vector *yt,            /* the current value y(t) of the sought function */
	double h,              /* the step to be taken */
	gsl_vector *yh,             /* output: y(t+h) */
	gsl_vector *err             /* output: error estimate */
);
int driver(
	void (*f)(double, gsl_vector*, gsl_vector*), /* right-hand-side of dy/dt=f(t,y) */
	double a,                     /* the start-point a */
	double b,                     /* the end-point of the integration */
	double h,                     /* initial step-size */
	double acc,                   /* absolute accuracy goal */
	double eps,                    /* relative accuracy goal */
	gsl_matrix *ylist,          //Matrix of stored y values 
    gsl_vector *xlist
); /* calculate y(b) */
#endif