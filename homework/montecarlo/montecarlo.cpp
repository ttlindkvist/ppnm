#include"montecarlo.hpp"
#include<stdlib.h>
#include<cmath>
#include<vector>
#include<assert.h>

std::pair<double,double> plainmc(double f(std::vector<double> &x), std::vector<double> &a, std::vector<double> &b, int N){
    assert(a.size() == b.size());
    int dim = a.size();
    double V=1;
    for(int i=0; i<dim; i++){ // Calculate volume of n-dimensional box
        V *= b[i] - a[i];
    }
    
    //Find point to sample within box
    std::vector<double> x = a;
    double sum = 0, sum2 = 0;
    for(int i=0; i<N; i++){
        for(int i=0; i<dim; i++){
            x[i] = 1.0*rand()/RAND_MAX*(b[i]-a[i]);
        }
        double fx = f(x);
        sum += fx;
        sum2 += fx*fx;
    }
    
    double mean = sum/N, sigma = sqrt(sum2/N-mean*mean);

    return std::pair<double, double>(mean*V, sigma*V/sqrt(N));
}
