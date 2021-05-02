#include<cstdio>
#include<cmath>
#include<gsl/gsl_matrix.h>
#include"minimization.hpp"
#include"function.hpp"

void print_vector(gsl_vector *v){
    int n = v->size;
    for(int i = 0; i<n; i++){
        printf("%10.7f ", gsl_vector_get(v, i));
    }
    printf("\n");
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

void partA(double eps){
    Function F(f1, 1);
    
    gsl_vector *x = gsl_vector_calloc(1);
    gsl_vector_set(x, 0, 2);
    
    printf("Minimizing x*x + 2x - 1\n");
    printf("Initial x   :   "); print_vector(x);
    int steps = qnewton(F, x, eps);
    printf("Minimum     :   "); print_vector(x);
    printf("Steps taken :    %d\n\n", steps);
    
    gsl_vector_free(x);
    
    Function rosenF(rosenbrock, 2);
    x = gsl_vector_alloc(2);
    gsl_vector_set(x, 0, 18);
    gsl_vector_set(x, 1, 15);
    
    printf("Minimizing Rosenbrock's function\n");
    printf("Initial x   :   "); print_vector(x);
    steps = qnewton(rosenF, x, eps);
    printf("Minimum     :   "); print_vector(x);
    printf("Steps taken :    %d\n\n", steps);
    
    
    Function himmelblauF(himmel, 2);
    gsl_vector_set(x, 0, 18);
    gsl_vector_set(x, 1, 15);
    
    printf("Minimizing Himmelblau's function\n");
    printf("Initial x   :   "); print_vector(x);
    steps = qnewton(himmelblauF, x, eps);
    printf("Minimum     :   "); print_vector(x);
    printf("Steps taken :    %d\n\n", steps);
    
    gsl_vector_free(x);
}

int main(){
    double eps = 1e-3;
    printf("Minimization\n--------- Part A ---------\nQuasi-newton method with symmetric Broyden update\neps=%g\n\n",eps);
    partA(eps);
    
    return 0;
}