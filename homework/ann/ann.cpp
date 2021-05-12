#include "ann.hpp"
#include "minimization.hpp"
#include<cmath>
#include<functional>
#include<cstdio>
#include<cstdlib>

ANN::ANN(int n, double (*f)(double)){
    this->params = gsl_vector_calloc(3*n);
    unsigned int seed = 1;
    for(int i = 0; i<3*n; i++){
        gsl_vector_set(this->params, i, 1.0*rand_r(&seed)/RAND_MAX);
    }
    this->activation_f = f; 
}

double ANN::cost(gsl_vector *p){
    gsl_vector_memcpy(this->params, p);
    double sum = 0;
    for(int i = 0; i<this->xs->size; i++){
        double x = gsl_vector_get(this->xs, i);
        sum += fabs(this->response(x) - gsl_vector_get(this->ys, i));
    }
    // printf("cost: %g\n", sum);
    return sum/xs->size;
}
void ANN::train(gsl_vector *xs, gsl_vector *ys){
    this->xs = xs;
    this->ys = ys;
    auto cost_lambda = [&](gsl_vector *p)->double {return this->cost(p);};
    std::function<double(gsl_vector*)> f = cost_lambda;
    gsl_vector *p = gsl_vector_alloc(this->params->size);
    gsl_vector_memcpy(p, this->params);
    int steps = qnewton(f, p, 1e-3);
    gsl_vector_memcpy(this->params, p);
    gsl_vector_free(p);
    
    printf("Training minimizing steps: %d\n", steps);
}
double ANN::response(double x){
    double sum = 0;
    for(int i = 0; i<this->params->size/3; i++){ //f((x-a)/b)*w
        sum += this->activation_f((x - gsl_vector_get(this->params, 3*i))/gsl_vector_get(this->params, 3*i+1))*gsl_vector_get(this->params, 3*i+2);
    }
    return sum;
}

ANN::~ANN(){
    gsl_vector_free(this->params);
}