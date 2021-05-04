#include "minimization.hpp"
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<cmath>
#include<cstdio>
#include<functional>

void print_vector(gsl_vector *v);

template<typename Callable> 
void gradient(Callable f, gsl_vector *grad, gsl_vector *x){
    static const double dx = sqrt(__DBL_EPSILON__);
    int n = x->size;
    
    double fx = f(x);
    for(int i = 0; i<n; i++){
        double xi = gsl_vector_get(x, i);
        
        gsl_vector_set(x, i, xi+dx);
        double f_new_x = f(x);
        gsl_vector_set(grad, i, (f_new_x-fx)/dx);
        
        gsl_vector_set(x, i, xi);
    }
}

template<typename Callable>
int qnewton(Callable f, gsl_vector *x, double eps){
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
    gradient(f, grad, x);
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
            gradient(f, y, x_s);
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
        gradient(f, grad, x);
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
template int qnewton<double(*)(gsl_vector*)>(double(*)(gsl_vector*), gsl_vector *, double);


double CurveFit::chi2(gsl_vector *params){
    double sum = 0;
    int size = this->xs.size();
    for(int i = 0; i<size; i++){
        double f = this->fit_function(this->xs[i], params);
        sum += (f-this->ys[i])*(f-this->ys[i])/this->yerrs[i]/this->yerrs[i];
    }
    return sum;
}

CurveFit::CurveFit(double (*f)(double, gsl_vector*), std::vector<double> &x, std::vector<double> &y, std::vector<double> &yerr)
    :xs(x),ys(y),yerrs(yerr){
   this->fit_function = f;
}
int CurveFit::fit(gsl_vector *params, double eps){
    int steps = qnewton([&](gsl_vector*p){return this->chi2(p);}, params, eps);
    return steps;
}
CurveFit::~CurveFit(){
}