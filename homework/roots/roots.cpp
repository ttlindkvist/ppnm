#include"roots.hpp"
#include"matrix.h"
#include"function.hpp"
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<cmath>
#include<assert.h>
#include<cstdio>

void newton(Function &f, gsl_vector *x, double eps, int &j_count, bool gradient){
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
        if(!gradient)   f.jacobian(J, x, fx);
        else            {f.hessian(J, x, fx);}

        j_count++;

        //Solve J \Delta x = -f(x) for \Delta x
        gsl_vector_scale(fx, -1);
        GS_decomp(J, R);
        GS_solve(J, R, fx, deltax);

        lambda = 2;
        do{
            lambda /= 2.;
            //Find norms
            fx_norm = gsl_blas_dnrm2(fx);

            gsl_vector_memcpy(x_ldx, x);
            gsl_blas_daxpy(lambda, deltax, x_ldx);

            if(!gradient)   f(x_ldx, fx_ldx);
            else            f.gradient(x_ldx, fx_ldx);

            fx_ldx_norm = gsl_blas_dnrm2(fx_ldx);
        } while(fx_ldx_norm > (1.-lambda/2)*fx_norm && lambda > 1./64);

        gsl_blas_daxpy(lambda, deltax, x);
    } while(fx_norm > eps);

    gsl_vector_free(fx);
    gsl_vector_free(fx_ldx);
    gsl_vector_free(deltax);
    gsl_vector_free(x_ldx);
    gsl_matrix_free(J);
    gsl_matrix_free(R);
}