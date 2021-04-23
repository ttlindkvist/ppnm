#ifndef __ROOTS_HPP__
#define __ROOTS_HPP__
#include<gsl/gsl_vector.h>

void newton(void f(gsl_vector *,gsl_vector *), gsl_vector *x, double eps, int &j_count);

#endif