#include<stdio.h>
#include<math.h>
#include"integration.h"

int fevals = 0; //Reused for all functions
double sqrt_x(double x){
    fevals++;
    return sqrt(x);
}
double circle(double x){
    fevals++;
    return 4.*sqrt(1-x*x);
}

int main(){
    fevals = 0;
    double I1 = adapt_quad24(sqrt_x, 0, 1, 1e-5, 1e-5);
    printf("int sqrt(x) from 0 to 1\nCalculated to \t%.10f\nShould be \t%.10f\nWith error \t%.10f\nwith %d function evaluations\n\n", \
            I1, 2./3, fabs(2./3-I1), fevals);
    
    fevals = 0;
    double I2 = adapt_quad24(circle, 0, 1, __DBL_EPSILON__, __DBL_EPSILON__);
    printf("int 4*sqrt(1-x*x) from 0 to 1\nCalculated to \t%.20f\nShould be \t%.20f\nWith error \t%.20f\nwith %d function evaluations\n\n", \
            I2, M_PI, fabs(M_PI-I2), fevals);
    fevals = 0;
    double I3 = adapt_clenshaw_curtis(sqrt_x, 0, 1, 1e-5, 1e-5);
    printf("int sqrt(x) from 0 to 1\nCalculated to \t%.10f\nShould be \t%.10f\nWith error \t%.10f\nwith %d function evaluations\n\n", \
            I3, 2./3, fabs(2./3-I3), fevals);
    
    fevals = 0;
    double I4 = adapt_clenshaw_curtis(circle, 0, 1, __DBL_EPSILON__, __DBL_EPSILON__);
    printf("int 4*sqrt(1-x*x) from 0 to 1\nCalculated to \t%.20f\nShould be \t%.20f\nWith error \t%.20f\nwith %d function evaluations\n\n", \
            I4, M_PI, fabs(M_PI-I4), fevals);
    return 0;
}