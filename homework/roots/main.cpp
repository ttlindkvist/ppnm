#include<cstdio>
#include<cmath>
#include<assert.h>
#include<functional>
#include<gsl/gsl_blas.h>
#include"roots.hpp"
#include"function.hpp"
#include"ode.h"


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

// -(1/2)f'' - (1/r)f = E*f
// f'' = -2*E*f - (2/r)f
static double _schr_E = 0; 
void schrodinger(double r, gsl_vector *y, gsl_vector *dydt){
    double f = gsl_vector_get(y, 0);
    double f_p = gsl_vector_get(y, 1);
    gsl_vector_set(dydt, 0, f_p);
    gsl_vector_set(dydt, 1, -2.0*(_schr_E + 1./r)*f);
}

static double _hydrogen_max_r = 8;
void hydrogen(gsl_vector *x, gsl_vector *fx){
    static double rmin = 1e-3;
    
    _schr_E = gsl_vector_get(x, 0);
    
    gsl_vector *fs = gsl_vector_alloc(2);

    //Boundary condition f(r->0) = r-r*r
    gsl_vector_set(fs, 0, rmin-rmin*rmin);
    gsl_vector_set(fs, 1, 1-2.*rmin);

    ode_driver_no_save(schrodinger, rmin, _hydrogen_max_r, 1e-3, 1e-3, 1e-3, 0.1, fs);

    gsl_vector_set(fx, 0, gsl_vector_get(fs, 0));

    gsl_vector_free(fs);
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

    printf("\n\n-------- PART B -------");

    double h_params[1] = {-2};
    gsl_vector_view h_view = gsl_vector_view_array(h_params, 1);

    Function hydrogen_energy(hydrogen, 1);
    printf("\n\nBound state of hydrogen via shooting method and initial guess of E=%g\n", h_params[0]);
    printf("eps = %g\n\n", eps);
    
    newton(hydrogen_energy, &h_view.vector, eps, j_count);

    printf("E = %.15g\nWith %d jacobi evaluations\n", h_params[0], j_count);
    printf("Actual solution\nE = -0.5\n");
}