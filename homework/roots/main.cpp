#include<cstdio>
#include<cmath>
#include<assert.h>
#include"roots.hpp"

void f1(gsl_vector *x, gsl_vector *fx){
    assert(x->size == fx->size && x->size ==2);
    double x1 = gsl_vector_get(x, 0);
    double x2 = gsl_vector_get(x, 1);
    gsl_vector_set(fx, 0, cos(x1)*sin(x2));
    gsl_vector_set(fx, 1, cos(x2));
}

int main(){
    double eps = 1e-4;
    int j_count = 0;
    double x[2] = {1, 1};
    
    gsl_vector_view x_view = gsl_vector_view_array(x, 2);

    printf("\nRoot finding\n");
    printf("-------- PART A -------\n\n");
    printf("Solving {cos(x)sin(y), cos(y)}=0 with initial guess (x,y) = (1, 1)\n");
    printf("eps = %g\n\n", eps);
    newton(f1, &x_view.vector, eps, j_count);
    printf("x1 = %.15g\nx2 = %.15g\nWith %d jacobi evaluations\n", x[0], x[1], j_count);
    printf("Actual solution\nx1 = %.15g\nx2 = %.15g\n", M_PI/2, M_PI/2);
}