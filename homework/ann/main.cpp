#include "ann.hpp"
#include "integration.h"
#include<cstdio>
#include<cmath>
#include<ccomplex>
#include<vector>

double gauss_wavelet(double x){
    return x*exp(-x*x);
}
double gauss_wavelet_deriv(double x){
    return exp(-x*x)*(1-2*x*x);
}
double gauss_wavelet_integ(double x){
    return -exp(-x*x)/2.;
}
double reLU(double x){
    return std::max(0.0, x);
}

double f(double x, void *p = NULL){
    return cos(5*x-1)*exp(-x*x);
}
double f_deriv(double x){
    return exp(-x*x)*(-5.*sin(5*x-1)-2*x*cos(5*x-1));
}
double f_integ(double a, double b){
    double err = 0;
    return adapt_quad24(f, a, b, 1e-3, 1e-3, &err);
}
void partAB(int n_neurons){
    std::vector<double> xs;
    std::vector<double> ys;
    double a = -1;
    double b = 1;
    int points = 20;
    for(int i = 0; i<points; i++){
        double x = a + (b-a)*i/(points-1);
        double y = f(x);
        xs.push_back(x);
        ys.push_back(y);
    }
    
    printf("------ PART A --------\n");
    printf("ANN with one hidden layer with %d neurons\n", n_neurons);
    printf("Done in C++ with an ANN class\n");
    printf("Activation function is a gaussian wavelet\n");
    printf("ANN is trained on a single set of points to fit the function cos(5x-1)*exp(-x*x)\n\n");
    
    ANN ann(n_neurons, gauss_wavelet, gauss_wavelet_deriv, gauss_wavelet_integ, xs, ys);
    int steps = ann.train();
    
    printf("Training took %d minimizing steps (with 10000 being max allowed pr. run)\n", steps);
    
    printf("\n------ PART B --------\n");
    printf("A plot of the derivative together with the integral over -1 to x has been - plotB.png\n");
    
    
    FILE *out = fopen("plot.out", "w");
    fprintf(out, "#index 0 -- points\n");
    for(int i = 0; i<xs.size(); i++){
        fprintf(out, "%g %g\n", xs[i], ys[i]);
    }
    
    fprintf(out, "\n\n#index 1 -- interpolated\n");
    int N = 100;
    for(int i = 0; i<=N; i++){
        double x = xs[0] + (xs.back() - xs[0])*1.0*i/N;
        fprintf(out, "%g %g\n", x, ann.response(x));
    }
    
    fprintf(out, "\n\n#index 2 -- derivatives and integrals\n");
    for(int i = 0; i<=N; i++){
        double x = xs[0] + (xs.back() - xs[0])*1.0*i/N;
        fprintf(out, "%g %g %g %g %g\n", x, ann.resp_deriv(x), ann.resp_integ(a, x), f_deriv(x), f_integ(a, x));
    }
    fclose(out);
}
int main(){
    printf("Artificial Neural Network Homework\n");
    partAB(4);
    
    return 0;
}