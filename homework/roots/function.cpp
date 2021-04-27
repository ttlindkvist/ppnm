#include"function.hpp"
#include<gsl/gsl_blas.h>
#include<cmath>
#include<assert.h>

Function::Function(double f(gsl_vector *), int n){
    f1 = f;
    f2 = NULL;
    this->n = n;
    _x_grad = gsl_vector_alloc(n);
    _x_jaco = gsl_vector_alloc(n);
    _f = gsl_vector_alloc(n);
}
Function::Function(void f(gsl_vector *, gsl_vector *), int n){
    f2 = f;
    f1 = NULL;
    this->n = n;
    _x_grad = gsl_vector_alloc(n);
    _x_jaco = gsl_vector_alloc(n);
    _f = gsl_vector_alloc(n);
}
void Function::gradient(gsl_vector *x, gsl_vector *grad){
    assert(f1 != NULL);
    static const double dx = sqrt(__DBL_EPSILON__);

    double fx = f1(x);
    for(int i = 0; i<n; i++){
        gsl_vector_memcpy(_x_grad, x);
        gsl_vector_set(_x_grad, i, gsl_vector_get(x, i)+dx);

        double f_new_x = f1(_x_grad);
        gsl_vector_set(grad, i, (f_new_x-fx)/dx);
    }
}
void Function::jacobian(gsl_matrix *J, gsl_vector *x, gsl_vector *fx){
    static const double dx = sqrt(__DBL_EPSILON__);
    assert(f2 != NULL);
    assert(J->size1 == J->size2 && J->size1 == x->size);
    int n = J->size1;

    f2(x, fx);

    for(int k = 0; k<n; k++){
        gsl_vector_memcpy(_x_jaco, x);
        gsl_vector_set(_x_jaco, k, gsl_vector_get(x, k)+dx);

        f2(_x_jaco, _f);
        for(int i = 0; i<n; i++){
            gsl_matrix_set(J, i, k, (gsl_vector_get(_f, i)-gsl_vector_get(fx, i))/dx);
        }
    }
}
//TODO: assuming continuity of second derivatives - the Hessian must be symmetric
void Function::hessian(gsl_matrix *J, gsl_vector *x, gsl_vector *fx){
    static const double dx = sqrt(__DBL_EPSILON__);
    assert(J->size1 == J->size2 && J->size1 == x->size);
    int n = J->size1;

    gradient(x, fx);
    
    for(int k = 0; k<n; k++){
        gsl_vector_memcpy(_x_jaco, x);
        gsl_vector_set(_x_jaco, k, gsl_vector_get(x, k)+dx);

        gradient(_x_jaco, _f);

        for(int i = 0; i<n; i++){
            gsl_matrix_set(J, i, k, (gsl_vector_get(_f, i)-gsl_vector_get(fx, i))/dx);
        }
    }
}
double Function::operator()(gsl_vector *x){
    assert(f1 != NULL);
    return f1(x);
}
void Function::operator()(gsl_vector *x, gsl_vector *fx){
    assert(f2 != NULL);
    f2(x, fx);
}
Function::~Function(){
    gsl_vector_free(_x_jaco);
    gsl_vector_free(_x_grad);
    gsl_vector_free(_f);
}