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
const int N = 6;
double c45[6] = {0, 1./4, 3./8, 12./13, 1, 1./2};
double b45[6] = {16./135, 0, 6656./12825, 28561./56430, -9./50, 2./55};
double b45e[6] = {25./216, 0, 1408./2565, 2197./4104, -1./5, 0};

// const int N = 2;
// double a12[1][1] =  {{1}};
// double c12[2] = {0, 1};
// double b12[2] = {0.5, 0.5};
// double b12e[2] = {1, 0};

void rkstep45(
	void (*f)(double, gsl_vector*, gsl_vector*), /* the f from dy/dt=f(t,y) */
	double t,                   /* the current value of the variable */
	gsl_vector *y,             /* the current value y(t) of the sought function */
	double h,                   /* the step to be taken */
	gsl_vector *yh,             /* output: y(t+h) */
	gsl_vector *err,             /* output: error estimate */
    gsl_matrix *Ks
){
    for(int i = 0; i<N; i++){
        gsl_vector_memcpy(yh, y);
        //Update predictor point
        for(int j = 0; j<i; j++){
            gsl_vector_view Kj = gsl_matrix_row(Ks, j);
            
            gsl_blas_daxpy(a45[i-1][j], &Kj.vector, yh);
        }
        //Calculate K_i
        gsl_vector_view Ki = gsl_matrix_row(Ks, i);
        f(t+c45[i]*h, yh, &Ki.vector);
        gsl_vector_scale(&Ki.vector, h);
        
        //Update error estimate
        gsl_blas_daxpy(b45[i]-b45e[i], &Ki.vector, err);
        
    }
    //Finally make y_n+1
    gsl_vector_memcpy(yh, y);
    for(int i = 0; i<N; i++){
        gsl_vector_view Ki = gsl_matrix_row(Ks, i);
        gsl_blas_daxpy(b45[i], &Ki.vector, yh);
    }
}
//For debugging
void rkstep12(
	void f(double t,gsl_vector*y,gsl_vector*dydt),
	double t,
	gsl_vector* y,
	double h,
	gsl_vector* yh,
	gsl_vector* err,
    gsl_matrix *K
){
	gsl_vector_view k0_view = gsl_matrix_row(K, 0);
	gsl_vector_view k1_view = gsl_matrix_row(K, 1);
    gsl_vector *k0 = &k0_view.vector;
    gsl_vector *k1 = &k1_view.vector;
    
	f(t,y,k0);
	gsl_vector_memcpy(yh,y);
	gsl_blas_daxpy(h/2,k0,yh);
	f(t+h/2,yh,k1);
	gsl_vector_memcpy(yh,y);
	gsl_blas_daxpy(h,k1,yh);
	gsl_vector_memcpy(err,k1);
	gsl_blas_daxpy(-1,k0,err);
	gsl_vector_scale(err,h);
}

int driver(
	void (*f)(double,gsl_vector*,gsl_vector*), /* right-hand-side of dy/dt=f(t,y) */
	double a,                     /* the start-point a */
	double b,                     /* the end-point of the integration */
	double h,                     /* initial step-size */
    double maxh,
	double acc,                   /* absolute accuracy goal */
	double eps,                   /* relative accuracy goal */
    gsl_matrix *ylist,            //Matrix of stored y values 
    gsl_vector *xlist
){
    int step = 0; //Current index of evaluation
    gsl_vector *err = gsl_vector_alloc(ylist->size2);
    gsl_matrix *Ks = gsl_matrix_calloc(N, ylist->size2);
    gsl_vector_view yt;
    gsl_vector_view yb;
    gsl_vector_set(xlist, 0, a);
    double x = a;
    while(x < b){
        gsl_vector_set_all(err, 0);
        yt = gsl_matrix_row(ylist, step);
        yb = gsl_matrix_row(ylist, step+1);
        if(x + h > b){
            h = b - x;
        }
        
        rkstep45(f, x, &yt.vector, h, &yb.vector, err, Ks);
        
        //Check if step is accepted
        double norme = gsl_blas_dnrm2(err); //Error at step 
        double normy = gsl_blas_dnrm2(&yb.vector); //Norm of yh
        double tol = (normy*eps+acc)*sqrt(h/(b-a)); //Tolerance this step
        if(norme<tol){ // Step is accepted
            step++;
            if(step > ylist->size1-2){ //If k grows larger than ylist matrix, make larger matrix
                return -step; //For now return -step
            }
            gsl_vector_set(xlist, step, x+h);
            x = x+h;
        }
        
        //Update step-size
        if(norme > 0){
            h *= pow(tol/norme, 0.25)*0.95;
        } else {h*=2;};
        if(h>maxh) h=maxh;
    }
    gsl_vector_free(err);
    gsl_matrix_free(Ks);
    return step+1;
}
int debug_driver(
	void (*f)(double,gsl_vector*,gsl_vector*), /* right-hand-side of dy/dt=f(t,y) */
	double a,                     /* the start-point a */
	double b,                     /* the end-point of the integration */
	double h,                     /* initial step-size */
    double maxh,
	double acc,                   /* absolute accuracy goal */
	double eps,                   /* relative accuracy goal */
    gsl_matrix *ylist,            //Matrix of stored y values 
    gsl_vector *xlist
){
    int step = 0; //Current index of evaluation
    gsl_vector *err = gsl_vector_alloc(ylist->size2);
    gsl_matrix *Ks = gsl_matrix_calloc(N, ylist->size2);
    gsl_vector_view yt;
    gsl_vector_view yb;
    gsl_vector_set(xlist, 0, a);
    double x = a;
    while(x < b){
        gsl_vector_set_all(err, 0);
        yt = gsl_matrix_row(ylist, step);
        yb = gsl_matrix_row(ylist, step+1);
        if(x + h > b){
            h = b - x;
        }
        fprintf(stderr, "test1\n");
        
        rkstep45(f, x, &yt.vector, h, &yb.vector, err, Ks);
        
        //Check if step is accepted
        double norme = gsl_blas_dnrm2(err); //Error at step 
        double normy = gsl_blas_dnrm2(&yb.vector); //Norm of yh
        double tol = (normy*eps+acc)*sqrt(h/(b-a)); //Tolerance this step
        fprintf(stderr, "test2\n");
        
        if(norme<tol){ // Step is accepted
            step++;
            fprintf(stderr, "test3 %d\n", step);
            if(step > ylist->size1-2){ //If k grows larger than ylist matrix, make larger matrix
                return -step; //For now return -step
            }
            
            gsl_vector_set(xlist, step, x+h);
            fprintf(stderr, "test12 %d\n", step);
            
            x = x+h;
        }
        fprintf(stderr, "test4\n");
        //Update step-size
        if(norme > 0){
            h *= pow(tol/norme, 0.25)*0.95;
        } else {h*=2;};
        if(h>maxh) h=maxh;
    }
    gsl_vector_free(err);
    gsl_matrix_free(Ks);
    return step+1;
}