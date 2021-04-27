#ifndef __ROOTS_HPP__
#define __ROOTS_HPP__
#include<gsl/gsl_vector.h>
#include"function.hpp"

void newton(Function &f, gsl_vector *x, double eps, int &j_count, bool gradient=false);

#endif