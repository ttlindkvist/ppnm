#ifndef __ANN_HPP__
#define __ANN_HPP__
#include<gsl/gsl_matrix.h>
#include<vector>

class ANN{
private:
    gsl_vector *params;
    int n_neurons; 
    std::vector<double> &xs;
    std::vector<double> &ys;
    double (*activation_f)(double);
    double (*activation_f_deriv)(double);
    double (*activation_f_integ)(double);
    double cost(gsl_vector *p);
public:
    ANN(int n, double(*f)(double), double(*f_deriv)(double), double(*f_integ)(double), std::vector<double> &xs, std::vector<double> &ys);
    int train();
    //Possibility of training with multiple datasets
    int train(std::vector<double> &xs, std::vector<double> &ys);
    double response(double x);
    double resp_deriv(double x);
    double resp_integ(double xa, double xb);
    ~ANN();
};

#endif