#ifndef __MINIMIZATION_HPP__
#define __MINIMIZATION_HPP__
#include<gsl/gsl_matrix.h>
#include<vector>

template<typename Callable>
extern int qnewton(Callable f, gsl_vector *x, double eps);

class CurveFit{
private:
    std::vector<double> &xs;
    std::vector<double> &ys;
    std::vector<double> &yerrs;
    double (*fit_function)(double, gsl_vector *);
public:
    CurveFit(double (*f)(double, gsl_vector*), std::vector<double> &x, std::vector<double> &y, std::vector<double> &yerr);
    double chi2(gsl_vector *params);
    void setPoints(std::vector<double> &x, std::vector<double> &y, std::vector<double> &yerr);
    int fit(gsl_vector *params, double eps);
    ~CurveFit();
};
    
    
#endif