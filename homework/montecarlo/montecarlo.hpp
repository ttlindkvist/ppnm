#ifndef __MONTECARLO_HPP__
#define __MONTECARLO_HPP__
#include<utility>
#include<vector>

std::pair<double,double> plainmc(double f(double *x), const std::vector<double> &a, const std::vector<double> &b, int N);
std::pair<double,double> quasi_mc(double f(double *x), const std::vector<double> &a, const std::vector<double> &b, int N);

#endif