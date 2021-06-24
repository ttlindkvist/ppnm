#include"matrix.h"
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>

#include<stdio.h>
#include<time.h>
#include<assert.h>
#include<complex.h>
#include<string.h>
#define true 1
#define false 0

//Testing algorithm on tall matrix with n>m
complex test_SVD(int n, int m, double tau, double eps){
    gsl_matrix *A = gsl_matrix_alloc(n, m);
    gsl_matrix *A_copy = gsl_matrix_alloc(n, m);
    gsl_matrix *D = gsl_matrix_alloc(m, m);
    gsl_matrix *U = gsl_matrix_alloc(n, m);
    gsl_matrix *V = gsl_matrix_alloc(m, m);
    gsl_matrix *UTU = gsl_matrix_alloc(m, m);
    gsl_matrix *VTV = gsl_matrix_alloc(m, m);
    
    gsl_matrix *P = gsl_matrix_alloc(n, m);
    gsl_matrix *tempM = gsl_matrix_alloc(n, m);
    
    
    gen_rand_matrix(A);
    gsl_matrix_memcpy(A_copy, A);
    int sweeps = SVD_two_jaco(A, D, V, U);
    
    //Check correct factorization
    gsl_matrix_memcpy(P, U);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, U, D, 0, tempM);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, tempM, V, 0, P);
    
    //Check orthogonality of U and V
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, U, U, 0, UTU);
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, V, V, 0, VTV);
    
    int err_code = 0;
    
    if(check_identity(UTU, tau, eps) == false){
        err_code = 1; // U is not orthogonal
    }
    else if(check_identity(VTV, tau, eps) == false){
        err_code = 2; // V is not orthogonal
    }
    else if(matrix_equals(P, A_copy, tau, eps) == false){
        err_code = 3; // Factorization is not equal to A
    }
    else if(check_diagonal(D, tau, eps) == false){
        err_code = 4; // Check D is diagonal
    }
        
    gsl_matrix_free(A);
    gsl_matrix_free(D);
    gsl_matrix_free(U);
    gsl_matrix_free(V);
    gsl_matrix_free(UTU);
    gsl_matrix_free(VTV);
    gsl_matrix_free(P);
    gsl_matrix_free(tempM);
    return err_code+I*sweeps;
}
void compare_GSL(){
    
}

void show_basic_func(){   
    printf("Two-sided Jacobi alg. for SVD\n\n");
    
    int n = 7, m = 5;
    gsl_matrix *A = gsl_matrix_alloc(n, m);    
    gsl_matrix *D = gsl_matrix_alloc(m, m);    
    gsl_matrix *U = gsl_matrix_alloc(n, m);    
    gsl_matrix *V = gsl_matrix_alloc(m, m);
    
    gsl_matrix *UTU = gsl_matrix_alloc(m, m);
    gsl_matrix *VTV = gsl_matrix_alloc(m, m);
    
    
    gsl_matrix *P = gsl_matrix_alloc(n, m);
    gsl_matrix *tempM = gsl_matrix_alloc(n, m);
    
    gen_rand_matrix(A);
    printf("Starting with %dx%d matrix: \n", n, m);
    print_matrix(A);
    
    int sweeps = SVD_two_jaco(A, D, V, U);
    printf("SVD done with %d sweeps\n", sweeps);
    
    printf("Check UDV^T=A: \n");
    gsl_matrix_memcpy(P, U);
    
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, U, D, 0, tempM);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, tempM, V, 0, P);
    
    print_matrix(P);
    
    printf("Matrix D: \n");
    print_matrix(D);
    printf("Matrix V: \n");
    print_matrix(V);
    printf("Matrix U: \n");
    print_matrix(U);
    
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, U, U, 0, UTU);
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, V, V, 0, VTV);
        
    printf("Checking orthonormality\nMatrix VT*V: \n");
    print_matrix(VTV);
    printf("Matrix UT*V: \n");
    print_matrix(UTU);
        
    gsl_matrix_free(A);
    gsl_matrix_free(D);
    gsl_matrix_free(U);
    gsl_matrix_free(V);
    gsl_matrix_free(UTU);
    gsl_matrix_free(VTV);
    gsl_matrix_free(P);
    gsl_matrix_free(tempM);    
}
void testing(){
    int n_tests = 10;
    int a = 50, b = 250;
    printf("\n\n----- TESTING -----\n");
    printf("Testing SVD algorithm on multiple random tall matrices, A, size NxM (N>M).\nM in [%d, %d] and N = M + rand([%d, %d]) + 1\n", a, b, a, b);
    printf("Testing is done by checking orthonormality of U and V, that D is diagonal, and U*D*V^T = A\n");
    for(int i = 0; i<n_tests; i++){
        int m = rand() % (b-a) + a;
        int n = m + rand() % (b-a) + 1;
        complex ret = test_SVD(n, m, 1e-5, 1e-5);
        int err_code = creal(ret);
        int sweeps = cimag(ret);
        printf("Test %d (size of A is %d x %d) with %d sweeps:  \t%s (code %d)\n", i+1, n, m, sweeps, err_code==0 ? "Success" : "Failure", err_code);
        
        if(err_code != 0){
            printf("ERROR with code %d\n", err_code);
            break;
        }
    }
}
void timing(){
    FILE *timing_output = fopen("timing.out", "w");
    
    const int runs = 15;
    const int max_N = 300;
    for(int n = max_N/runs; n<=max_N; n+=max_N/runs){
        fprintf(stderr, "timing n=%d\n", n);
        gsl_matrix *A = gsl_matrix_alloc(n, n);    
        gsl_matrix *Acopy = gsl_matrix_alloc(n, n);
        gsl_vector *S = gsl_vector_alloc(n);
        gsl_matrix *U = gsl_matrix_alloc(n, n);    
        gsl_matrix *V = gsl_matrix_alloc(n, n);
        gen_rand_matrix(A);
        gsl_matrix_memcpy(Acopy, A);
        
        clock_t start = clock();
        SVD_two_jaco_square(A, V, U);
        clock_t end = clock();
        int msec_my_algo = (double)(end-start) * 1000. / CLOCKS_PER_SEC;
        
        start = clock();
        gsl_linalg_SV_decomp_jacobi(Acopy, V, S);
        end = clock();
        int msec_gsl = (double)(end-start) * 1000. / CLOCKS_PER_SEC;
        
        
        fprintf(timing_output, "%d %d %d\n", n, msec_my_algo, msec_gsl);
        
        gsl_matrix_free(A);
        gsl_matrix_free(Acopy);
        gsl_matrix_free(U);
        gsl_matrix_free(V);
        gsl_vector_free(S);
    }
    fclose(timing_output);
}

int main(int argc, char**argv){
    srand(42);
    
    printf("EXAMINATION PROJECT 13\n");
    printf("AUTHOR: Thomas Toft Lindkvist\n");    
    printf("AUID643642 - Student number: 201905635\n");
    printf("35 mod 22 = 13\n"); 
    
    if(argc >= 2){
        if(strcmp(argv[1], "test") == 0){
            show_basic_func();
            testing();
        } else if(strcmp(argv[1], "time") == 0){
            timing();
        } else {
            fprintf(stderr, "Unknown commandline argument passed\n");
            return EXIT_FAILURE;
        }
    } else{
        show_basic_func();
    }    
    return EXIT_SUCCESS; 
}