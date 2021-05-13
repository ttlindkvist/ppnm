#include "ann.hpp"
#include "minimization.hpp"
#include<cmath>
#include<functional>
#include<cstdio>
#include<cstdlib>

ANN::ANN(int n, double (*f)(double), double(*f_deriv)(double), double(*f_integ)(double), std::vector<double> &xs, std::vector<double> &ys) : xs(xs), ys(ys){
    this->n_neurons = n;
    this->activation_f = f; 
    this->activation_f_deriv = f_deriv; 
    this->activation_f_integ = f_integ; 
    
    this->params = gsl_vector_calloc(3*n);
    
    unsigned int seed = 2;
    double a = xs.front();
    double b = xs.back();
    for(int i = 0; i<n; i++){
        //Place offsets uniformly over interval
        gsl_vector_set(this->params, 3*i, a + (b-a)*i/(n-1) + (b-a)*(2.0*rand_r(&seed)/RAND_MAX-1)*1e-2);
        gsl_vector_set(this->params, 3*i+1, 1.0);
        gsl_vector_set(this->params, 3*i+2, 1.0);
    }
}

double ANN::cost(gsl_vector *p){
    gsl_vector_memcpy(this->params, p);
    double n = this->xs.size();
    
    double sum = 0;
    for(int i = 0; i<n; i++){;
        sum += fabs(this->response(this->xs[i]) - this->ys[i]);
    }
    // printf("cost: %g\n", sum);
    return sum/this->xs.size();
}
int ANN::train(){
    //Lambda function with capture of object by reference
    auto cost_lambda = [&](gsl_vector *p)->double {return this->cost(p);};
    //Recast lambda to c++ functional, necessary for use of qnewton template function
    std::function<double(gsl_vector*)> f = cost_lambda; 
    
    gsl_vector *p = gsl_vector_alloc(this->params->size);
    gsl_vector_memcpy(p, this->params);
    
    int steps = qnewton(f, p, 1e-3);
    
    gsl_vector_memcpy(this->params, p);
    gsl_vector_free(p);
    
    return steps;
}
int ANN::train(std::vector<double> &xs, std::vector<double> &ys){
    this->xs = xs;
    this->ys = ys;
    return this->train();    
}
double ANN::response(double x){
    double sum = 0;
    double a = 0;
    double b = 0;
    double w = 0;
    for(int i = 0; i<this->n_neurons; i++){ //f((x-a)/b)*w
        a = gsl_vector_get(this->params, 3*i + 0);
        b = gsl_vector_get(this->params, 3*i + 1);
        w = gsl_vector_get(this->params, 3*i + 2);
        sum += this->activation_f((x - a)/b)*w;
    }
    return sum;
}
double ANN::resp_deriv(double x){
    double sum = 0;
    double a = 0;
    double b = 0;
    double w = 0;
    for(int i = 0; i<this->n_neurons; i++){
        a = gsl_vector_get(this->params, 3*i + 0);
        b = gsl_vector_get(this->params, 3*i + 1);
        w = gsl_vector_get(this->params, 3*i + 2);
        
        sum += this->activation_f_deriv((x - a)/b)*w/b;
    }    
    return sum;
}
double ANN::resp_integ(double xa, double xb){
    double sum = 0;
    double a = 0;
    double b = 0;
    double w = 0;
    for(int i = 0; i<this->n_neurons; i++){
        a = gsl_vector_get(this->params, 3*i + 0);
        b = gsl_vector_get(this->params, 3*i + 1);
        w = gsl_vector_get(this->params, 3*i + 2);
        
        sum += this->activation_f_integ((xa - a)/b)*w*b - this->activation_f_integ((xb - a)/b)*w*b;
    }    
    return sum;
}


ANN::~ANN(){
    gsl_vector_free(this->params);
}