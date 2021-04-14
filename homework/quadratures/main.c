#include<stdio.h>
#include<math.h>
#include"integration.h"
#include <gsl/gsl_integration.h>

int fevals = 0; //Reused for all functions
double sqrt_x(double x, void *params){
    fevals++;
    return sqrt(x);
}
double one_over_sqrt(double x, void *params){
        fevals++;
        return 1./sqrt(x);
}
double ln_over_sqrt(double x, void *params){
    fevals++;
    return log(x)/sqrt(x);
}
double circle(double x, void *params){
    fevals++;
    return 4.*sqrt(1-x*x);
}
void partA(double eps, double abs){
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
}
void partB(double eps, double abs, double eps2, double abs2){
    printf("\n------ PART B ------\n\n");
    
    fevals = 0;
    double I = adapt_clenshaw_curtis(one_over_sqrt, 0, 1, abs, eps);
    printf("int 1/sqrt(x) from 0 to 1 - using Clenshaw-Curtis transformation\nCalculated to \t%.16f\nShould be \t%.16f\nWith error \t%.15f\nwith %d function evaluations\n\n", \
            I, 2., fabs(2.-I), fevals);
    
    fevals = 0;
    I = adapt_clenshaw_curtis(ln_over_sqrt, 0, 1, abs, eps);
    printf("int ln(x)/sqrt(x) from 0 to 1 - using Clenshaw-Curtis transformation\nCalculated to \t%.16f\nShould be \t%.16f\nWith error \t %.15f\nwith %d function evaluations\n\n", \
            I, -4., fabs(-4.-I), fevals);
    
    printf("As seen the CC transformation results in much fewer function evaluations\n");
    printf("Also, the error seems to be on the same order or lower (at a first glance)\n");

    printf("\nEvaluating 4*sqrt(1-x*x) with and without Clenshaw-Curtis with\n");
    printf("Absolute precision %g - Relative precision %g\n", abs2, eps2);
    fevals = 0;
    I = adapt_quad24(circle, 0, 1, abs2, eps2);
    printf("int 4*sqrt(1-x*x) from 0 to 1\nCalculated to \t%.16f\nShould be \t%.16f\nWith error \t%.16f\nwith %d function evaluations\n\n", \
            I, M_PI, fabs(M_PI-I), fevals);

    fevals = 0;
    I = adapt_clenshaw_curtis(circle, 0, 1, abs2, eps2);
    printf("int 4*sqrt(1-x*x) from 0 to 1 - using CC transformation\nCalculated to \t%.16f\nShould be \t%.16f\nWith error \t%.16f\nwith %d function evaluations\n\n", \
            I, M_PI, fabs(M_PI-I), fevals);   
    printf("As seen the CC transformation results in lower overall error, even with the same abs and eps supplied, with around the same number of function evaluations\n");
}
void partB_GSL(double eps, double abs){
    printf("\n\nComparing to the GSL implementation\n\n");
    gsl_integration_workspace *w = gsl_integration_workspace_alloc (1000);
    double result, error;

    gsl_function F;
    F.function = &ln_over_sqrt;

    fevals = 0;
    gsl_integration_qags(&F, 0, 1, abs, eps, 1000,
                            w, &result, &error);
    printf("int ln(x)/sqrt(x) from 0 to 1\n");
    printf ("Calculated to   %.16f\n", result);
    printf ("exact result    %.16f\n", -4.);
    printf ("estimated error  %.16f\n", error);
    printf ("actual error     %.16f\n", fabs(-4.-result));
    printf ("intervals        %zu\n", w->size);
    printf ("calls            %i\n", fevals);

  gsl_integration_workspace_free (w);
}
void partC(){

}

int main(){
    double eps = 1e-4, abs = 1e-4;
    partA(eps, abs);
    partB(eps, abs, 1e-10, 1e-10);
    partB_GSL(eps, abs);
    partC();

    return 0;
}