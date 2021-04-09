#include"integration.h"
#include<assert.h>
#include<math.h>

/*
  Points used are {1/6, 2/6, 4/6, 5/6}
  Weights {2/6, 1/6, 1/6, 2/6} (trapezium rule)
  Weights for embedded error-estimate {1/4, 1/4, 1/4, 1/4} (rectangle rule)
  
*/

double n_points = 4;
double absc24[] = {1./6, 2./6, 4./6, 5./6};
double w24[] = {2./6, 1./6, 1./6, 2./6}; 
double v24[] = {1./4, 1./4, 1./4, 1./4}; 

double qag24(double f(double), double a, double b, double abs, double eps, double f2, double f3, unsigned int rec_depth){
    assert(rec_depth<1000000); // Make sure the depth of recursion is below some max depth
    
    double f1 = f(a + (b-a)*absc24[0]);
    double f4 = f(a + (b-a)*absc24[3]);
    
    //Evaluate 2nd and 4th order estimate
    double Q = (2*f1 + f2 + f3 + 2*f4)/6*(b-a);
    double q = (f1 + f2 + f3 + f4)/4*(b-a);
    
    double tolerance = abs + eps*fabs(Q);
    double error = fabs(Q-q);

    if(error < tolerance){ // Accept evaluation
        return Q;
    } else { //Subdivide interval
        double midpoint = (a+b)/2;
        double Q1 = qag24(f, a, midpoint, abs/sqrt(2.), eps, f1, f2, rec_depth + 1);
        double Q2 = qag24(f, midpoint, b, abs/sqrt(2.), eps, f3, f4, rec_depth + 1);
        return Q1+Q2;
    }
    
}
double adapt_quad24(double f(double), double a, double b, double abs, double eps){
    //Calculate the middle two points
    double f2 = f(a + (b-a)*absc24[1]);
    double f3 = f(a + (b-a)*absc24[2]);
    return qag24(f, a, b, abs, eps, f2, f3, 0);    
}