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

double gamma_integrand(double x, void *params){
    double z = *(double*)params;
    return pow(x, z-1.)*exp(-x);
}
double integ_gamma(double z){
    if(z<0) return integ_gamma(z+1)/z;
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
    double result, error;
    gsl_function F;
    F.function = &gamma_integrand;
    F.params = &z;
    gsl_integration_qagiu(&F, 0, 1e-3, 1e-5, 1000, w, &result, &error);
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
    FILE *gamma_out = fopen("gamma.dat", "w");

    double xmin=-4, xmax = 4;
    int N = 1000;
    for(int i = 0; i<N; i++){
        double x = xmin + (xmax-xmin)/N*i;
        fprintf(erf_out, "%.10f \t %.10f \t %.10f\n", x, integ_erf(x), gsl_sf_erf(x));
    }

    xmin=-5, xmax = 7;
    for(int i = 0; i<N; i++){
        double x = xmin + (xmax-xmin)/N*i + 0.003;
        fprintf(gamma_out, "%.10f \t %.10f \t %.10f\n", x, integ_gamma(x), gsl_sf_gamma(x));
    }

    fclose(erf_out);
    fclose(gamma_out);
    return 0;
}