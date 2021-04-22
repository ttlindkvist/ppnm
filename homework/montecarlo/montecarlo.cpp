#include"montecarlo.hpp"
#include<stdlib.h>
#include<cmath>
#include<assert.h>
#include<omp.h>
#include<iostream>

std::pair<double,double> plainmc(double f(double *x), const std::vector<double> &a, const std::vector<double> &b, int N){
    assert(a.size() == b.size());
    
    int dim = a.size();
    double V = 1;

    double diff[dim];
    for(int i=0; i<dim; i++){ // Calculate volume of n-dimensional box
        diff[i] = b[i] - a[i];
        V *= diff[i];
    }
    
    //Find point to sample within box
    double x[dim];
    double sum = 0, sum2 = 0;
    for(int i = 0; i<N; i++){
        for(int j = 0; j<dim; j++){
            x[j] = a[j] + 1.0*rand()/RAND_MAX*diff[j];
        }
        double fx = f(x);
        sum += fx;
        sum2 += fx*fx;
    }
    
    double mean = sum/N, sigma = sqrt(sum2/N-mean*mean);

    return std::pair<double, double>(mean*V, sigma*V/sqrt(N));
}

//Quasi-random generators
double corput(int n, int base){
    double q = 0, bk = 1./base;
    while(n > 0){
        q += (n % base)*bk;
        n /= base;
        bk /= base;
    }
    return q;
}
void halton1(int n, int d, double *x){
    static const int base[]= {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67};
    static const int maxd = sizeof(base)/sizeof(int);
    
    assert(d<=maxd);
    for(int i = 0; i<d; i++){
        x[i] = corput(n, base[i]);
    }
}
void halton2(int n, int d, double *x){
    static const int base[]= {3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71};
    static const int maxd = sizeof(base)/sizeof(int);
    
    assert(d<=maxd);
    for(int i = 0; i<d; i++){
        x[i] = corput(n, base[i]);
    }
}

std::pair<double,double> quasi_mc(double f(double *x), const std::vector<double> &a, const std::vector<double> &b, int N){
    assert(a.size() == b.size());
    
    int dim = a.size();
    double V = 1;

    double diff[dim];
    for(int i=0; i<dim; i++){ // Calculate volume of n-dimensional box
        diff[i] = b[i] - a[i];
        V *= diff[i];
    }

    double x1[dim];
    double x2[dim];
    
    double sum1=0, sum2=0;
    #pragma omp parallel sections
    {
        #pragma omp section //halton1
        {
            for(int i = 0; i<N/2; i++){
                halton1(i+1, dim, x1);
                for(int j = 0; j<dim; j++){
                    x1[j] = x1[j]*diff[j]+a[j];
                }
                sum1 += f(x1);
            }
        }
        #pragma omp section //halton2
        {
            for(int i = 0; i<N/2; i++){
                halton2(i+1, dim, x2);
                for(int j = 0; j<dim; j++){
                    x2[j] = x2[j]*diff[j]+a[j];
                }
                sum2 += f(x2);
            }
        }
    }
    return std::pair<double,double>((sum1+sum2)/N*V, fabs(sum1-sum2)/N*V);
}