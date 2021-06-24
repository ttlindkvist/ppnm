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

static double _hydrogen_min_r = 1e-3;
static double _hydrogen_max_r = 8;
void hydrogen(gsl_vector *x, gsl_vector *fx){
    _schr_E = gsl_vector_get(x, 0);
    
    gsl_vector *fs = gsl_vector_alloc(2);

    //Boundary condition f(r->0) = r-r*r
    gsl_vector_set(fs, 0, _hydrogen_min_r-_hydrogen_min_r*_hydrogen_min_r);
    gsl_vector_set(fs, 1, 1-2.*_hydrogen_min_r);

    ode_driver_no_save(schrodinger, _hydrogen_min_r, _hydrogen_max_r, 1e-3, 1e-3, 1e-3, 0.1, fs);

    gsl_vector_set(fx, 0, gsl_vector_get(fs, 0));

    gsl_vector_free(fs);
}
void hydrogen_improved_boundary(gsl_vector *x, gsl_vector *fx){
    _schr_E = gsl_vector_get(x, 0);
    
    gsl_vector *fs = gsl_vector_alloc(2);

    //Boundary condition f(r->0) = r-r*r
    gsl_vector_set(fs, 0, _hydrogen_min_r-_hydrogen_min_r*_hydrogen_min_r);
    gsl_vector_set(fs, 1, 1-2.*_hydrogen_min_r);

    ode_driver_no_save(schrodinger, _hydrogen_min_r, _hydrogen_max_r, 1e-3, 1e-3, 1e-3, 0.1, fs);

    //Improved boundary condition for f(r->infty) = r*e^(-kr), k=sqrt(-2E) gives correction to final function value
    gsl_vector_set(fx, 0, gsl_vector_get(fs, 0) - _hydrogen_max_r*exp(-sqrt(-2*_schr_E)*_hydrogen_max_r));

    gsl_vector_free(fs);
}
void partA(double eps){
    int j_count = 0;
    double x[2] = {1, 1};
    
    gsl_vector_view x_view = gsl_vector_view_array(x, 2);

    Function f(f1, 2);
    
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
void partB(double eps){
    _hydrogen_min_r = 1e-3;
    _hydrogen_max_r = 8;
    
    int j_count = 0;
    double h_params[1] = {-2};
    gsl_vector_view h_view = gsl_vector_view_array(h_params, 1);

    Function hydrogen_energy(hydrogen, 1);
    printf("\n\nBound state of hydrogen via shooting method and initial guess of E=%g\n", h_params[0]);
    printf("eps = %g\n\n", eps);
    
    newton(hydrogen_energy, &h_view.vector, eps, j_count);

    printf("E = %.15g\nWith %d jacobi evaluations\n", h_params[0], j_count);
    printf("Actual solution\nE = -0.5\n");

    dyn_matrix *ylist = dyn_matrix_alloc(20, 2);
    dyn_vector *xlist = dyn_vector_alloc(20);
    dyn_matrix_set(ylist, 0, 0, _hydrogen_min_r - _hydrogen_min_r*_hydrogen_min_r);
    dyn_matrix_set(ylist, 0, 1, 1 - 2*_hydrogen_min_r);
    dyn_vector_set(xlist, 0, _hydrogen_min_r);

    int steps = ode_driver_save(schrodinger, _hydrogen_min_r, _hydrogen_max_r, 1e-3, 1e-4, 1e-4, 0.01, ylist, xlist);

    FILE *output = fopen("wavefunc.out", "w");

    for(int i = 0; i<steps; i++){
        double x = dyn_vector_get(xlist, i);
        double *y = dyn_matrix_row(ylist, i);
        double exact = x*exp(-x);
        fprintf(output, "%g %g %g %g\n", x, y[0], y[1], exact);
    }

    fclose(output);
    dyn_matrix_free(ylist);
    dyn_vector_free(xlist);

    printf("\nThe calculated and exact wavefunction have been plotted to wavefunc.png\n");
}
void partC(double eps){
    Function hydrogen_energy(hydrogen, 1);  
    Function hydrogen_energy_improved(hydrogen_improved_boundary, 1);
    int steps = 100;
    double r_max_start = 1;
    double r_max_end = 10;
    
    for(int i = 0; i<steps; i++){
        double r_max = r_max_start + (r_max_end-r_max_start)/steps*i;  
        
        _hydrogen_max_r = r_max;
        int j_count = 0;
        int j_count2 = 0;
        double h_params[1] = {-2};

        gsl_vector_view h_view = gsl_vector_view_array(h_params, 1);
        newton(hydrogen_energy, &h_view.vector, eps, j_count);
        double E_old_bound = h_params[0];
        
        h_params[0] = -2;        
        newton(hydrogen_energy_improved, &h_view.vector, eps, j_count2);
        double E_new_bound = h_params[0];
        
        printf("rmax = %g - E_old = %.15g - E_new = %.15g\n", r_max, E_old_bound, E_new_bound);
    }
}


int main(){
    printf("\nRoot finding\n");
    
    // printf("-------- PART A -------\n\n");
    // partA(1e-4);

    // printf("\n\n-------- PART B -------");
    // partB(1e-4);

    printf("\n\n-------- PART C -------\n");
    partC(1e-4);
}