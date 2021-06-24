#include"matrix.h"
#include<gsl/gsl_blas.h>
#include<gsl/gsl_eigen.h>

#include<stdio.h>
#include<assert.h>

int testDiag(gsl_matrix *A){
    return 0;
}

int testOrtho(gsl_matrix *A){
    return 0;
}


int main(){
    srand(42);
    printf("EXAMINATION PROJECT 13\n");
    printf("AUTHOR: Thomas Toft Lindkvist\n");    
    printf("AUID643642 - Student number: 201905635\n");
    printf("35 mod 22 = 13\n");    
    printf("Two-sided Jacobi alg. for SVD\n\n");
    
    
    int n = 5;
    gsl_matrix *A = gsl_matrix_alloc(n, n);    
    gsl_matrix *U = gsl_matrix_alloc(n, n);    
    gsl_matrix *V = gsl_matrix_alloc(n, n);    
    gsl_matrix *P = gsl_matrix_alloc(n, n);
    gsl_matrix *tempM = gsl_matrix_alloc(n, n);
    gen_rand_matrix(A);
    printf("Starting with: \n");
    print_matrix(A);
    
    int sweeps = SVD_two_jaco_square(A, V, U);
    printf("SVD done with %d sweeps\n", sweeps);
    
    printf("Check UDV^T=A: \n");
    gsl_matrix_memcpy(P, U);
    
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, U, A, 0, tempM);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, tempM, V, 0, P);
    
    print_matrix(P);
    
    printf("Matrix D: \n");
    print_matrix(A);
    printf("Matrix V: \n");
    print_matrix(V);
    printf("Matrix U: \n");
    print_matrix(U);
        
    
    gsl_matrix_free(A);
    gsl_matrix_free(U);
    gsl_matrix_free(V);
    gsl_matrix_free(P);
    gsl_matrix_free(tempM);
    
    
    printf("TODO:\n");
    printf("Testing on random NxN matrices with N=[100, 300]\n");
    printf("Testing is done by checking orthogonality of U and V, diagonality of D, and comparing eigenvalues with a gsl-implementation\n\n");
    return 0;    
}