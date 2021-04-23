#include"roots.hpp"
#include"matrix.h"
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<cmath>
#include<assert.h>
#include<cstdio>

void jacobian(void f(gsl_vector *x,gsl_vector *fx), gsl_matrix *J, gsl_vector *x, gsl_vector *fx){
    static const double dx = sqrt(__DBL_EPSILON__);
    assert(J->size1 == J->size2 && J->size1 == x->size);
    int n = J->size1;

    f(x, fx);

    gsl_vector *new_x = gsl_vector_alloc(n);
    gsl_vector *fi = gsl_vector_alloc(n);

    for(int k = 0; k<n; k++){
        gsl_vector_set_basis(new_x, k);
        gsl_vector_scale(new_x, dx);
        gsl_blas_daxpy(1., x, new_x);
        f(new_x, fi);
        for(int i = 0; i<n; i++){
            gsl_matrix_set(J, i, k, (gsl_vector_get(fi, i)-gsl_vector_get(fx, i))/dx);
        }
    }
    gsl_vector_free(fi);
    gsl_vector_free(new_x);
}


void newton(void f(gsl_vector *x,gsl_vector *fx), gsl_vector *x, double eps, int &j_count){
    gsl_vector *fx = gsl_vector_alloc(x->size);
    gsl_vector *fx_ldx = gsl_vector_alloc(x->size);
    gsl_vector *deltax = gsl_vector_alloc(x->size);
    gsl_vector *x_ldx = gsl_vector_alloc(x->size);
    gsl_matrix *J = gsl_matrix_alloc(x->size, x->size);
    gsl_matrix *R = gsl_matrix_calloc(x->size, x->size);

    double fx_norm = 0;
    double fx_ldx_norm = 0;
    double lambda = 1;
    j_count = 0;
    do{
        jacobian(f, J, x, fx);
        j_count++;

        //Solve J \Delta x = -f(x) for \Delta x
        gsl_vector_scale(fx, -1);
        GS_decomp(J, R);
        GS_solve(J, R, fx, deltax);

        lambda = 1;
        //Find norms
        fx_norm = gsl_blas_dnrm2(fx);

        gsl_vector_memcpy(x_ldx, x);
        gsl_blas_daxpy(lambda, deltax, x_ldx);
        f(x_ldx, fx_ldx);
        fx_ldx_norm = gsl_blas_dnrm2(fx_ldx);


        while(fx_ldx_norm > (1.-lambda/2)*fx_norm && lambda > 1./64){ //Perhaps switch for do-while
            lambda /= 2.;
            //Find norms
            fx_norm = gsl_blas_dnrm2(fx);

            gsl_vector_memcpy(x_ldx, x);
            gsl_blas_daxpy(lambda, deltax, x_ldx);
            f(x_ldx, fx_ldx);
            fx_ldx_norm = gsl_blas_dnrm2(fx_ldx);
        }
        gsl_blas_daxpy(lambda, deltax, x);

    } while(fx_norm > eps);

    gsl_vector_free(fx);
    gsl_vector_free(fx_ldx);
    gsl_vector_free(deltax);
    gsl_vector_free(x_ldx);
    gsl_matrix_free(J);
    gsl_matrix_free(R);
}
