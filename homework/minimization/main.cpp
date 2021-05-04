#include<cstdio>
#include<cmath>
#include<assert.h>
#include<vector>
#include<gsl/gsl_matrix.h>
#include"minimization.hpp"

void print_vector(gsl_vector *v){
    int n = v->size;
    for(int i = 0; i<n; i++){
        printf("%10.7f ", gsl_vector_get(v, i));
    }
    printf("\n");
}
void read_points(const char* filename, std::vector<double> &xs, std::vector<double> &ys, std::vector<double> &y_errs){
    FILE *in_stream = fopen(filename, "r");
    double x, y, yerr;
    while(fscanf(in_stream, "%lg %lg %lg", &x, &y, &yerr) > 0){
        xs.push_back(x);
        ys.push_back(y);
        y_errs.push_back(yerr);
    }
    fclose(in_stream);
}

double f1(gsl_vector *x){
    double x0 = gsl_vector_get(x, 0);
    return x0*x0 + 2*x0 - 1; 
}
double rosenbrock(gsl_vector *x_vec){
    double x = gsl_vector_get(x_vec, 0);
    double y = gsl_vector_get(x_vec, 1);
    return (1.-x)*(1.-x) + 100.*(y-x*x)*(y-x*x);
}
double himmel(gsl_vector *x_vec){
    double x = gsl_vector_get(x_vec, 0);
    double y = gsl_vector_get(x_vec, 1);
    return pow((x*x+y-11), 2) + (x+y*y-7)*(x+y*y-7);
}
double BW(double E, gsl_vector *params){
    double m = gsl_vector_get(params, 0);
    double Gamm = gsl_vector_get(params, 1);
    double A = gsl_vector_get(params, 2);
    return A/((E-m)*(E-m) + Gamm*Gamm/4);   
}

void partA(double eps){    
    printf("Minimization\n--------- Part A ---------\nQuasi-newton method with symmetric Broyden update\neps=%g\n\n",eps);
    gsl_vector *x = gsl_vector_calloc(1);
    gsl_vector_set(x, 0, 2);
    
    printf("Minimizing x*x + 2x - 1\n");
    printf("Initial x   :   "); print_vector(x);
    int steps = qnewton(f1, x, eps);
    printf("Minimum     :   "); print_vector(x);
    printf("Steps taken :    %d\n\n", steps);
    
    gsl_vector_free(x);
    
    x = gsl_vector_alloc(2);
    gsl_vector_set(x, 0, 18);
    gsl_vector_set(x, 1, 15);
    
    printf("Minimizing Rosenbrock's function\n");
    printf("Initial x   :   "); print_vector(x);
    steps = qnewton(rosenbrock, x, eps);
    printf("Minimum     :   "); print_vector(x);
    printf("Steps taken :    %d\n\n", steps);
    
    
    gsl_vector_set(x, 0, 18);
    gsl_vector_set(x, 1, 15);
    
    printf("Minimizing Himmelblau's function\n");
    printf("Initial x   :   "); print_vector(x);
    steps = qnewton(himmel, x, eps);
    printf("Minimum     :   "); print_vector(x);
    printf("Steps taken :    %d\n\n", steps);
    
    gsl_vector_free(x);
}
void partB(double eps){
    printf("\n\n----------- Part B -----------\neps=%g\n\n",eps);
    
    std::vector<double> xs, ys, yerrs;
    read_points("higgs.in", xs, ys, yerrs);
    
    gsl_vector *p = gsl_vector_calloc(3);
    gsl_vector_set(p, 0, 120);
    gsl_vector_set(p, 1, 4);
    gsl_vector_set(p, 2, 5);
    
    printf("Fitting higgs data to BW function\nInitial guess for (m, Gamma, A) = ");
    print_vector(p);
    
    CurveFit fit(BW, xs, ys, yerrs);
    int steps = fit.fit(p, 1e-3);
    
    printf("\nFitting yields (m, Gamma, A) = ");    
    print_vector(p);
    printf("with %d steps taken\n", steps);
    
    printf("Note: Gamma has a wrong sign, since BW has Gamma*Gamma - positive is of course used\n");
    
    FILE *higgs_out = fopen("higgs_fit.out", "w");
    double E1=100, E2=160;
    int N = 1000;
    for(int i = 0; i<N; i++){
        double E = (E2-E1)*i/N + E1;
        double f = BW(E, p);
        fprintf(higgs_out, "%g %g\n", E, f);
    }
    fclose(higgs_out);
    
    gsl_vector_free(p);
}
int main(){
    double eps = 1e-3;
    partA(eps);
    partB(eps);
    
    return 0;
}