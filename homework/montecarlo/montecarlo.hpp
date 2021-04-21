#ifndef __MONTECARLO_HPP__
#define __MONTECARLO_HPP__
#include<utility>
#include<vector>

std::pair<double,double> plainmc(double f(std::vector<double> &x), std::vector<double> &a, std::vector<double> &b, int N);

#endif