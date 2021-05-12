#ifndef __ANN_HPP__
#define __ANN_HPP__
#include<gsl/gsl_matrix.h>

class ANN{
private:
    gsl_vector *params;
    gsl_vector *xs;
    gsl_vector *ys;
    double (*activation_f)(double);
    double cost(gsl_vector *p);
public:
    ANN(int n, double(*f)(double));
    void train(gsl_vector *xs, gsl_vector *ys);
    double response(double x);
    ~ANN();
};

#endif