#include<stdio.h>
#include<math.h>
#include<gsl/gsl_integration.h>
#include<gsl/gsl_sf_erf.h>
#include<gsl/gsl_sf_gamma.h>

double f(double x, void *params){
    return log(x)/sqrt(x);
}
double erf_integrand(double t, void *params){
    return 2.0/M_SQRTPI*exp(-t*t);
}
double integ_erf(double z){
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
    double result, error;
    gsl_function F;
    F.function = &erf_integrand;
    gsl_integration_qags(&F, 0, z, 0, 1e-7, 1000, w, &result, &error);
    gsl_integration_workspace_free(w);

    return result;
}

typedef struct {double x; double n;} bessel_params;

double J_integrand(double tau, void *p){
    bessel_params *params = (bessel_params*)p;
    double x = params->x;
    double n = params->n;

    return cos(n*tau-x*sin(tau))/M_PI;
}
double integ_J(double x, int n){
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
    double result, error;
    gsl_function F;
    F.function = &J_integrand;
    bessel_params params = {x, n};
    F.params = &params;

    gsl_integration_qags(&F, 0, M_PI, 1e-9, 1e-6, 1000, w, &result, &error);
    gsl_integration_workspace_free(w);
    
    return result;
}

int main(){
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
    double result, error;
    gsl_function F;
    F.function = &f;
    gsl_integration_qags(&F, 0, 1, 0, 1e-7, 1000, w, &result, &error);
    printf("Part A: result = %.10f\terror = %.10g\n", result, error);
    gsl_integration_workspace_free(w);

    FILE *erf_out = fopen("erf.dat", "w");
    FILE *bessel_out = fopen("bessel.dat", "w");

    double xmin=-4, xmax = 4;
    int N = 1000;
    for(int i = 0; i<N; i++){
        double x = xmin + (xmax-xmin)/N*i;
        fprintf(erf_out, "%.10f \t %.10f \t %.10f\n", x, integ_erf(x), gsl_sf_erf(x));
    }

    xmin=0, xmax = 20;
    for(int i = 0; i<N; i++){
        double x = xmin + (xmax-xmin)/N*i;
        fprintf(bessel_out, "%.10f \t %.10f \t %.10f \t %.10f\n", x, integ_J(x, 0), integ_J(x, 1), integ_J(x, 2));
    }

    fclose(erf_out);
    fclose(bessel_out);
    return 0;
}