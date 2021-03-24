#include"matrix.h"
#include<stdio.h>


int main(){
    int n = 4;
    gsl_matrix *A = gsl_matrix_alloc(n, n);
    gsl_matrix *V = gsl_matrix_alloc(n, n);
    gen_rand_symm_matrix(A);
    
    print_matrix(A);
    
    jacobi_diag(A, V);
    
    print_matrix(A);
    print_matrix(V);
    
    
    gsl_matrix_free(A);
    gsl_matrix_free(V);
    return 0;
}