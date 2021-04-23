#include"integration.h"
#include<assert.h>
#include<math.h>
#define NULL ((void*)0)
/*
  Points used are {1/6, 2/6, 4/6, 5/6}
  Weights {2/6, 1/6, 1/6, 2/6} (trapezium rule)
  Weights for embedded error-estimate {1/4, 1/4, 1/4, 1/4} (rectangle rule)
  
*/

double n_points = 4;
double absc24[] = {1./6, 2./6, 4./6, 5./6};
double w24[] = {2./6, 1./6, 1./6, 2./6}; 
double v24[] = {1./4, 1./4, 1./4, 1./4}; 

double qag24(double f(double,void*), double a, double b, double abs, double eps, double f2, double f3, unsigned int rec_depth, double *err){
    assert(rec_depth<1000000); // Make sure the depth of recursion is below some max depth
    

    double f1 = f(a + (b-a)*absc24[0], NULL);
    double f4 = f(a + (b-a)*absc24[3], NULL);
    
    //Evaluate 2nd and 4th order estimate
    double Q = (2*f1 + f2 + f3 + 2*f4)/6*(b-a);
    double q = (f1 + f2 + f3 + f4)/4*(b-a);
    
    double tolerance = abs + eps*fabs(Q);
    double error = fabs(Q-q);

    if(error < tolerance){ // Accept evaluation
        *err = error;
        return Q;
    } else { //Subdivide interval
        double err1, err2;
        double midpoint = (a+b)/2;
        double Q1 = qag24(f, a, midpoint, abs/sqrt(2.), eps, f1, f2, rec_depth + 1, &err1);
        double Q2 = qag24(f, midpoint, b, abs/sqrt(2.), eps, f3, f4, rec_depth + 1, &err2);
        *err = sqrt(err1*err1 + err2*err2); // Calculate combined error from sub-intervals
        return Q1+Q2;
    }
    
}
double adapt_quad24(double f(double,void*), double a, double b, double abs, double eps, double *err){
    //Calculate the middle two points
    double f2 = f(a + (b-a)*absc24[1], NULL);
    double f3 = f(a + (b-a)*absc24[2], NULL);
    //Reset error estimate
    double I = qag24(f, a, b, abs, eps, f2, f3, 0, err);
    return I;
}

static double (*cc_func)(double, void*);
static double cc_a,cc_b;
double cc_wrapper(double x, void* params){
    return cc_func((cc_a+cc_b)/2 + (cc_b-cc_a)*cos(x)/2, NULL)*(cc_b-cc_a)*sin(x)/2;
}

double adapt_clenshaw_curtis(double f(double,void*), double a, double b, double abs, double eps, double *err){
    cc_func = f;
    cc_a = a;
    cc_b = b;
    
    return adapt_quad24(cc_wrapper, 0, M_PI, abs, eps, err);  
}

static double (*inf_func)(double, void*);
static double inf_a,inf_b;
double inf_inf_wrapper(double x, void* params){ //From -1 to 1
    return inf_func(x/(1-x*x), NULL)*(1+x*x)/((1-x*x)*(1-x*x));
}
double a_inf_wrapper(double x, void* params){ // From 0 to 1
    return inf_func(inf_a + (1-x)/x, NULL)/(x*x);
}
double inf_b_wrapper(double x, void* params){ // From 0 to 1
    return inf_func(inf_b - (1-x)/x,NULL)/(x*x);
}

double adapt_inf(double f(double,void*), double a, double b, double abs, double eps, double *err){
    inf_func = f;
    inf_a = a;
    inf_b = b;
    
    if(isinf(a) && isinf(b)){
        return adapt_clenshaw_curtis(inf_inf_wrapper, -1, 1, abs, eps, err);
    }
    else if(isinf(b)){
        return adapt_clenshaw_curtis(a_inf_wrapper, 0, 1, abs, eps, err);
    }
    else if(isinf(a)){
        return adapt_clenshaw_curtis(inf_b_wrapper, 0, 1, abs, eps, err);
    }
    return adapt_clenshaw_curtis(f, a, b, abs, eps, err);
}