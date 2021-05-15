#ifndef __INTEGRATION_H__
#define __INTEGRATION_H__

double qag24(double f(double, void*), double a, double b, double abs, double eps, double f2, double f3, unsigned int rec_depth, double *err);
double adapt_quad24(double f(double, void*), double a, double b, double abs, double eps, double *err);
double adapt_clenshaw_curtis(double f(double, void*), double a, double b, double abs, double eps, double *err);
double adapt_inf(double f(double,void*), double a, double b, double abs, double eps, double *err);


#endif