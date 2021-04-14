#ifndef __INTEGRATION_H__
#define __INTEGRATION_H__

double qag24(double f(double), double a, double b, double abs, double eps, double f2, double f3, unsigned int rec_depth);
double adapt_quad24(double f(double), double a, double b, double abs, double eps);
double adapt_clenshaw_curtis(double f(double), double a, double b, double abs, double eps);

#endif