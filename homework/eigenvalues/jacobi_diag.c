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
void jacobi_diag(gsl_matrix *A, gsl_matrix *V){
    assert(A->size1 == A->size2);
    int n = A->size1;
    int changes = 0;
    
    //Initialize V to identity matrix
    gsl_matrix_set_identity(V);

    //Sweep over upper triangle of A untill rotation doesn't change anything
    do{
        for(int p = 0; p<n-1; p++){
            for(int q = p+1; q<n; q++){
                //Check diagonal elements
                double Apq = gsl_matrix_get(A, p, q);
                double App = gsl_matrix_get(A, p, p);
                double Aqq = gsl_matrix_get(A, q, q);
                double theta = 0.3*atan2(2*Apq, Aqq-App);
                double c = cos(theta), s = sin(theta);
                double new_App = c*c*App - 2*s*c*Apq + s*s*Aqq;
                double new_Aqq = s*s*App + 2*s*c*Apq + c*c*Aqq;
                if(equals(new_App, App, 1e-7, 1e-7)==0 || equals(new_Aqq, Aqq, 1e-7, 1e-7)==0){
                    changes++;
                    timesJ(A, p, q, theta);
                    Jtimes(A, p, q, -theta);
                    timesJ(V, p, q, theta);
                }
            }
        }
    } while(changes > 0);
}