#ifndef __FUNCTION_HPP__
#define __FUNCTION_HPP__
#include<gsl/gsl_matrix.h>

class Function{
    private:
        double (*f1)(gsl_vector*);
        void (*f2)(gsl_vector*,gsl_vector*);
        gsl_vector *_x_jaco;
        gsl_vector *_f;
        int n;
    public:
        Function(double f(gsl_vector *), int n);
        Function(void f(gsl_vector *,gsl_vector *), int n);
        void gradient(gsl_vector *grad, gsl_vector *x);
        void jacobian(gsl_matrix *J, gsl_vector *x, gsl_vector *fx);
        void hessian(gsl_matrix *H, gsl_vector *x, gsl_vector *fx);
        double operator()(gsl_vector *);
        void operator()(gsl_vector *, gsl_vector *);
        ~Function();
};

#endif