#include<cstdio>
#include<cmath>
#include<assert.h>
#include<functional>
#include<gsl/gsl_blas.h>
#include"roots.hpp"
#include"function.hpp"


void f1(gsl_vector *x, gsl_vector *fx){
    assert(x->size == fx->size && x->size ==2);
    double x1 = gsl_vector_get(x, 0);
    double x2 = gsl_vector_get(x, 1);
    gsl_vector_set(fx, 0, cos(x1)*sin(x2));
    gsl_vector_set(fx, 1, cos(x2));
}
double rosenbrock(gsl_vector *xs){
    double x = gsl_vector_get(xs, 0);
    double y = gsl_vector_get(xs, 1);
    return (1-x)*(1-x) + 100*(y-x*x)*(y-x*x);
}
void rosenbrock_gradient(gsl_vector *x, gsl_vector *grad){
    static const double dx = sqrt(__DBL_EPSILON__);

    gsl_vector *new_x = gsl_vector_alloc(2);

    double fx = rosenbrock(x);
    for(int i = 0; i<2; i++){
        gsl_vector_set_basis(new_x, i);
        gsl_vector_scale(new_x, dx);
        gsl_blas_daxpy(1., x, new_x);
        double f_new_x = rosenbrock(new_x);
        gsl_vector_set(grad, i, (f_new_x-fx)/dx);
    }
    gsl_vector_free(new_x);
}

int main(){
    double eps = 1e-4;
    int j_count = 0;
    double x[2] = {1, 1};
    
    gsl_vector_view x_view = gsl_vector_view_array(x, 2);

    Function f(f1, 2);
    printf("\nRoot finding\n");
    printf("-------- PART A -------\n\n");
    printf("Solving {cos(x)sin(y), cos(y)}=0 with initial guess (x,y) = (1, 1)\n");
    printf("eps = %g\n\n", eps);
    newton(f, &x_view.vector, eps, j_count);
    printf("x1 = %.15g\nx2 = %.15g\nWith %d jacobi evaluations\n", x[0], x[1], j_count);
    printf("Actual solution\nx1 = %.15g\nx2 = %.15g\n", M_PI/2, M_PI/2);


    x[0] = 0.5;
    x[1] = 0.5;
    Function rosenb(rosenbrock, 2);
    printf("\n\nFinding extrema of rosenbrock with initial guess (0.5, 0.5)\n");
    printf("eps = %g\n\n", eps);
    newton(rosenb, &x_view.vector, eps, j_count, true);
    printf("x1 = %.15g\nx2 = %.15g\nWith %d hessian evaluations\n", x[0], x[1], j_count);
    printf("Actual solution\nx1 = %.15g\nx2 = %.15g\n", 1., 1.);
}