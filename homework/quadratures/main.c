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
void test_function_24(double f(double, void*), double a, double b, double abs, double eps, double exact){
    fevals = 0;
    double error_estimate = 0;
    double I = adapt_quad24(f, 0, 1, abs, eps, &error_estimate);
    printf("Calculated to     %.16f\n", I);
    printf("Exact result      %.16f\n", exact);
    printf("Estimated error   %.16f\n", error_estimate);
    printf("Actual error      %.16f\n", fabs(exact-I));
    printf("with %d function evaluations\n\n", fevals);
}
void test_function_CC(double f(double, void*), double a, double b, double abs, double eps, double exact){
    fevals = 0;
    double error_estimate = 0;
    double I = adapt_clenshaw_curtis(f, 0, 1, abs, eps, &error_estimate);
    printf("Calculated to     %.16f\n", I);
    printf("Exact result      %.16f\n", exact);
    printf("Estimated error   %.16f\n", error_estimate);
    printf("Actual error      %.16f\n", fabs(exact-I));
    printf("with %d function evaluations\n\n", fevals);
}

void partA(double eps, double abs){
    printf("\n------ PART A ------ \nAll integrals calculated with recursive adaptive integrator\n");
    printf("Absolute precision %g - Relative precision %g\n", abs, eps);

    printf("int sqrt(x) from 0 to 1\n");
    test_function_24(sqrt_x, 0, 1, abs, eps, 2./3);

    printf("int 4*sqrt(1-x*x) from 0 to 1\n");
    test_function_24(circle, 0, 1, abs, eps, M_PI);

    printf("int 1/sqrt(x) from 0 to 1\n");
    test_function_24(one_over_sqrt, 0, 1, abs, eps, 2.);
    
    printf("int ln(x)/sqrt(x) from 0 to 1\n");
    test_function_24(ln_over_sqrt, 0, 1, abs, eps, -4.);
}
void partB(double eps, double abs, double eps2, double abs2){
    printf("\n------ PART B ------\n\n");
    
    printf("int 1/sqrt(x) from 0 to 1 - using Clenshaw-Curtis transformation\n");
    test_function_CC(one_over_sqrt, 0, 1, abs, eps, 2.);

    printf("int ln(x)/sqrt(x) from 0 to 1 - using Clenshaw-Curtis transformation\n");
    test_function_CC(ln_over_sqrt, 0, 1, abs, eps, -4.);

    printf("As seen the CC transformation results in much fewer function evaluations\n");
    printf("Also, the error seems to be on the same order or lower (at a first glance)\n");

    printf("\nEvaluating 4*sqrt(1-x*x) with and without Clenshaw-Curtis with\n");
    printf("Absolute precision %g - Relative precision %g\n\n", abs2, eps2);

    printf("int 4*sqrt(1-x*x) from 0 to 1\n");
    test_function_24(circle, 0, 1, abs2, eps2, M_PI);

    printf("int 4*sqrt(1-x*x) from 0 to 1 - using CC transformation\n");
    test_function_CC(circle, 0, 1, abs2, eps2, M_PI);

    printf("As seen the CC transformation results in lower overall error, even with the same abs and eps supplied, with around the same number of function evaluations\n");
}
void partB_GSL(double eps, double abs){
    printf("\n\nComparing to the GSL implementation\n\n");
    gsl_integration_workspace *w = gsl_integration_workspace_alloc (1000);
    double result, error;

    gsl_function F;

    F.function = &one_over_sqrt;
    fevals = 0;
    gsl_integration_qags(&F, 0, 1, abs, eps, 1000, w, &result, &error);
    printf("int 1/sqrt(x) from 0 to 1\n");
    printf("Calculated to    %.16f\n", result);
    printf("Exact result     %.16f\n", 2.);
    printf("Estimated error  %.16f\n", error);
    printf("Actual error     %.16f\n", fabs(2.-result));
    printf("with %i function evaluations\n\n", fevals);

    F.function = &ln_over_sqrt;
    fevals = 0;
    gsl_integration_qags(&F, 0, 1, abs, eps, 1000, w, &result, &error);
    printf("int ln(x)/sqrt(x) from 0 to 1\n");
    printf("Calculated to   %.16f\n", result);
    printf("Exact result    %.16f\n", -4.);
    printf("Estimated error  %.16f\n", error);
    printf("Actual error     %.16f\n", fabs(-4.-result));
    printf("with %i function evaluations\n\n", fevals);
    
    gsl_integration_workspace_free (w);

    printf("Even though the GSL implementation uses more integrand evaluations, the error is many orders of magnitude lower.\n");
    printf("Despite the abs and eps supplied are the same\n");
}
void partC(){
    printf("\n------ PART C ------\n\n");
    printf("As seen above the implemented integrator returns its estimate of the integration error\n");
    printf("On should also note that this estimate always will be larger than the actual error\n");
}

int main(){
    double eps = 1e-4, abs = 1e-4;
    partA(eps, abs);
    partB(eps, abs, 1e-10, 1e-10);
    partB_GSL(eps, abs);
    partC();

    return 0;
}