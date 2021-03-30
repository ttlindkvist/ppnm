#include"ode.h"
#include<gsl/gsl_blas.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<math.h>

double a45[5][5] = {
    {1./4,  0, 0, 0, 0},
    {3./32, 9./32, 0, 0, 0},
    {1932./2197, -7200./2197, 7296./2197, 0, 0},
    {439./216, -8, 3680./513, -845./4104, 0},
    {-8./27, 2, -3544./2565, 1859./4104, -11./40}
};
double c45[6] = {0, 1./4, 3./8, 12./13, 1, 1/2};
double b45[6] = {16./135, 0, 6656./12825, 28561./56430, -9./50, 2./55};
double b45e[6] = {25./216, 0, 1408./2565, 2197./4104, -1./5, 0};

void rkstep45(
	void (*f)(double, gsl_vector*, gsl_vector*), /* the f from dy/dt=f(t,y) */
	double t,                   /* the current value of the variable */
	gsl_vector* yt,             /* the current value y(t) of the sought function */
	double h,                   /* the step to be taken */
	gsl_vector* yh,             /* output: y(t+h) */
	gsl_vector* err             /* output: error estimate */
){
    gsl_vector *K = gsl_vector_alloc(yt->size);
    gsl_vector *yb[6];
    for(int i = 0; i<6; i++){
        yb[i] = gsl_vector_alloc(yt->size);
        gsl_vector_memcpy(yb[i], yt);
    }
    gsl_vector_memcpy(yh, yt);
    
    for(int i = 0; i<6; i++){
        //Calculate K_i
        f(t+c45[i]*h, yb[i], K);
        gsl_vector_scale(K, h);
        
        //Update y_{n+1} and error estimate
        gsl_blas_daxpy(b45[i], K, yh);
        gsl_blas_daxpy(b45[i]-b45e[i], K, err);
        
        //Update next proping pointgs yb
        for(int j = i+1; j<6; j++){
            gsl_blas_daxpy(a45[j][i], K, yb[j]);
        }
    }
    
    gsl_vector_free(K);
    for(int i = 0; i<6; i++){
        gsl_vector_free(yb[i]);
    }
}

int driver(
	void (*f)(double,gsl_vector*,gsl_vector*), /* right-hand-side of dy/dt=f(t,y) */
	double a,                     /* the start-point a */
	double b,                     /* the end-point of the integration */
	double h,                     /* initial step-size */
	double acc,                   /* absolute accuracy goal */
	double eps,                   /* relative accuracy goal */
    gsl_matrix *ylist,            //Matrix of stored y values 
    gsl_vector *xlist
){
    int k = 0; //Current index of evaluation
    double x;
    gsl_vector *err = gsl_vector_alloc(ylist->size2);
    gsl_vector_set(xlist, 0, a);
    while(gsl_vector_get(xlist, k)<b){
        gsl_vector_set_all(err, 0);
        x = gsl_vector_get(xlist, k);
        gsl_vector_view yt = gsl_matrix_row(ylist, k);
        gsl_vector_view yb = gsl_matrix_row(ylist, k+1);
        if(x+h > b){
            h = b-x;
        }
        rkstep45(f, x, &yt.vector, h, &yb.vector, err);
        
        double ei = gsl_blas_dnrm2(err); //Error at step 
        double normy = gsl_blas_dnrm2(&yb.vector); //Norm of yh
        double tol = (normy*eps+acc)*sqrt(h/(b-a)); //Tolerance this step
        if(ei<tol){ // Step is accepted
            k++;
            if(k>ylist->size1-2){ //If k grows larger than ylist matrix, make larger matrix
                return -k; //For now return -1
            }
            gsl_vector_set(xlist, k, x+h);
        }
        if(ei > 0){
            h *= pow(tol/ei, 0.25)*0.95;
        } else {h*=2;};
    }
    gsl_vector_free(err);
    return k+1;
}