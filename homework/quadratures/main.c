#include<stdio.h>
#include<math.h>
#include"integration.h"

int fevals = 0; //Reused for all functions
double sqrt_x(double x){
    fevals++;
    return sqrt(x);
}
double one_over_sqrt(double x){
        fevals++;
        return 1./sqrt(x);
}
double ln_over_sqrt(double x){
    fevals++;
    return log(x)/sqrt(x);
}
double circle(double x){
    fevals++;
    return 4.*sqrt(1-x*x);
}

int main(){
    double eps = 1e-4, abs = 1e-4;

    printf("\n------ PART A ------ \nAll integrals calculated with recursive adaptive integrator\n");
    printf("Absolute precision %g - Relative precision %g\n", abs, eps);

    fevals = 0;
    double I = adapt_quad24(sqrt_x, 0, 1, abs, eps);
    printf("int sqrt(x) from 0 to 1\nCalculated to \t%.16f\nShould be \t%.16f\nWith error\t%.16f\nwith %d function evaluations\n\n", \
            I, 2./3, fabs(2./3-I), fevals);

    fevals = 0;
    I = adapt_quad24(circle, 0, 1, abs, eps);
    printf("int 4*sqrt(1-x*x) from 0 to 1\nCalculated to \t%.16f\nShould be \t%.16f\nWith error\t%.16f\nwith %d function evaluations\n\n", \
            I, M_PI, fabs(M_PI-I), fevals);

    fevals = 0;
    I = adapt_quad24(one_over_sqrt, 0, 1, abs, eps);
    printf("int  1/sqrt(x) from 0 to 1\nCalculated to \t%.16f\nShould be \t%.16f\nWith error \t%.16f\nwith %d function evaluations\n\n", \
            I, 2., fabs(2.-I), fevals);
    
    fevals = 0;
    I = adapt_quad24(ln_over_sqrt, 0, 1, abs, eps);
    printf("int ln(x)/sqrt(x) from 0 to 1\nCalculated to \t%.16f\nShould be \t%.16f\nWith error \t %.16f\nwith %d function evaluations\n\n", \
            I, -4., fabs(-4.-I), fevals);
    
    printf("\n------ PART B ------\n\n");
    
    fevals = 0;
    I = adapt_clenshaw_curtis(one_over_sqrt, 0, 1, abs, eps);
    printf("int 1/sqrt(x) from 0 to 1 - using Clenshaw-Curtis transformation\nCalculated to \t%.16f\nShould be \t%.16f\nWith error \t%.15f\nwith %d function evaluations\n\n", \
            I, 2., fabs(2.-I), fevals);
    
    fevals = 0;
    I = adapt_clenshaw_curtis(ln_over_sqrt, 0, 1, abs, eps);
    printf("int ln(x)/sqrt(x) from 0 to 1 - using Clenshaw-Curtis transformation\nCalculated to \t%.16f\nShould be \t%.16f\nWith error \t %.15f\nwith %d function evaluations\n\n", \
            I, -4., fabs(-4.-I), fevals);
    
    printf("As seen the CC transformation results in much fewer function evaluations\n");
    printf("Also, the error seems to be on the same order or lower (at a first glance)\n");

    printf("\nEvaluating 4*sqrt(1-x*x) with and without Clenshaw-Curtis with machine epsilon as eps and abs tolerance\n");
    fevals = 0;
    I = adapt_quad24(circle, 0, 1, __DBL_EPSILON__, __DBL_EPSILON__);
    printf("int 4*sqrt(1-x*x) from 0 to 1\nCalculated to \t%.16f\nShould be \t%.16f\nWith error \t%.16f\nwith %d function evaluations\n\n", \
            I, M_PI, fabs(M_PI-I), fevals);

    fevals = 0;
    I = adapt_clenshaw_curtis(circle, 0, 1, __DBL_EPSILON__, __DBL_EPSILON__);
    printf("int 4*sqrt(1-x*x) from 0 to 1\nCalculated to \t%.16f\nShould be \t%.16f\nWith error \t%.16f\nwith %d function evaluations\n\n", \
            I, M_PI, fabs(M_PI-I), fevals);   

    return 0;
}