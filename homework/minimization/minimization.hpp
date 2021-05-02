#ifndef __MINIMIZATION_HPP__
#define __MINIMIZATION_HPP__
#include"function.hpp"
#include<gsl/gsl_matrix.h>
int qnewton(Function &f, gsl_vector *x, double eps);


#endif