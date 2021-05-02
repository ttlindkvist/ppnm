#include "minimization.hpp"
#include"function.hpp"
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<cmath>
#include<cstdio>


int qnewton(Function &f, gsl_vector *x, double eps){
    const double ALPHA = 1e-3;
    const double DELTA = sqrt(__DBL_EPSILON__); 
    int n = x->size;
    gsl_matrix *B = gsl_matrix_alloc(n, n); //Inverse Hessian
    gsl_matrix_set_identity(B);

    gsl_vector *s = gsl_vector_alloc(n);
    gsl_vector *grad = gsl_vector_alloc(n);
    
    gsl_vector *x_s = gsl_vector_alloc(n);
    gsl_vector *y = gsl_vector_alloc(n);
    gsl_vector *u = gsl_vector_alloc(n);
    gsl_vector *a = gsl_vector_alloc(n);

    double fx;
    double lamb = 1;
    
    int steps = 0;
    
    fx = f(x);
    f.gradient(grad, x);
        
    while(gsl_blas_dnrm2(grad) > eps){        
        gsl_blas_dgemv(CblasNoTrans, -1, B, grad, 0, s);
        
        //Find appropriate step-size
        lamb = 2;
        double fz = 0;
        double sTgrad = 0;
        do{
            //Update fz
            gsl_vector_memcpy(x_s, x);
            gsl_blas_daxpy(1, s, x_s); // Calc x_s = x + s
            fz = f(x_s);
            
            
            gsl_blas_ddot(s, grad, &sTgrad); 
            
            //Update s for next loop
            gsl_vector_scale(s, 0.5);
            lamb /= 2.;            
        }while(fz > fx + ALPHA*sTgrad && lamb >= DELTA);
        
        //Update B - Hessian inverse
        if(lamb < DELTA){
            gsl_matrix_set_identity(B);
        } else {
            //Symmetric Broyden update
            f.gradient(y, x_s);
            gsl_blas_daxpy(-1, grad, y);
            
            gsl_blas_dgemv(CblasNoTrans, -1, B, y, 0, u);
            gsl_blas_daxpy(1, s, u);
            
            double sTy; gsl_blas_ddot(s, y, &sTy);
            if(sTy > DELTA){
                double uTy; gsl_blas_ddot(u, y, &uTy);
                gsl_vector_memcpy(a, u);
                gsl_blas_daxpy(-uTy/2./sTy, s, a);
                gsl_vector_scale(a, 1./sTy);
                
                //Perform update
                gsl_blas_dger(1, a, s, B);
                gsl_blas_dger(1, s, a, B);
            } else {
                gsl_matrix_set_identity(B);
            }
        }
        //Update x
        gsl_vector_memcpy(x, x_s);
        steps++;
        fx = f(x);
        f.gradient(grad, x);
    }
    
    gsl_matrix_free(B);
    gsl_vector_free(grad);
    gsl_vector_free(s);
    gsl_vector_free(x_s);
    gsl_vector_free(y);
    gsl_vector_free(u);
    gsl_vector_free(a);
    return steps;
}