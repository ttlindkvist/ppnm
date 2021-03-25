#include"matrix.h"
#include<gsl/gsl_blas.h>
#include<math.h>
#include<assert.h>

void gen_rand_symm_matrix(gsl_matrix *A){
    assert(A->size1 == A->size2);
    int n = A->size1;
    for(int i = 0; i<n; i++){
        for(int j = i; j<n; j++){
            double r = 1.0*rand()/RAND_MAX;
            gsl_matrix_set(A, i, j, r);
            gsl_matrix_set(A, j, i, r);
        }
    }
}

void timesJ(gsl_matrix *A, int p, int q, double theta){
    double c = cos(theta), s = sin(theta);
    //Changes the columns p and q in A
    for(int i = 0; i<A->size1; i++){
        double Aip = gsl_matrix_get(A, i, p);
        double Aiq = gsl_matrix_get(A, i, q);
        gsl_matrix_set(A, i, p, c*Aip-s*Aiq);
        gsl_matrix_set(A, i, q, s*Aip+c*Aiq);
    }
}
void Jtimes(gsl_matrix *A, int p, int q, double theta){
    double c = cos(theta), s = sin(theta);
    //Changes the rows p and q in A
    for(int j = 0; j<A->size2; j++){
        double Apj = gsl_matrix_get(A, p, j);
        double Aqj = gsl_matrix_get(A, q, j);
        gsl_matrix_set(A, p, j,  c*Apj + s*Aqj);
        gsl_matrix_set(A, q, j, -s*Apj + c*Aqj);
    }
}
double sum_off_diagonal(gsl_matrix *A){
    int n = A->size1;
    double sum=0;
    for(int i = 0; i<n; i++){
        for(int j = i+1; j<n; j++){
            sum += fabs(gsl_matrix_get(A, i, j));
        }
    }
    return sum;
}
void jacobi_diag(gsl_matrix *A, gsl_matrix *V){
    assert(A->size1 == A->size2);
    int n = A->size1;
    // int changes = 0;
    
    //Initialize V to identity matrix
    gsl_matrix_set_identity(V);

    //Sweep until sum of off-diagonal elements are smaller than some eps=double machine epsilon
    do{
        for(int p = 0; p<n-1; p++){
            for(int q = p+1; q<n; q++){
                
                double Apq = gsl_matrix_get(A, p, q);
                double App = gsl_matrix_get(A, p, p);
                double Aqq = gsl_matrix_get(A, q, q);
                double theta = 0.5*atan2(2*Apq, Aqq-App);
                
                timesJ(A, p, q, theta);
                Jtimes(A, p, q, -theta);
                timesJ(V, p, q, theta);
            }
        }
    } while(sum_off_diagonal(A)>__DBL_EPSILON__);
}