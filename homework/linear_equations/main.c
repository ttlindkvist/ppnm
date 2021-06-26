#include"matrix.h"
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>
#include<stdio.h>
#include<assert.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

//Returns 0 on success
//Rerturns non-zero and different for each possible failure
int test_partA1(int n, int m, int print){
    int return_code = 0;
    
    gsl_matrix *A = gsl_matrix_alloc(n, m);
    gsl_matrix *R = gsl_matrix_calloc(m, m);
    gsl_matrix *Q = gsl_matrix_alloc(n, m);
    gsl_matrix *QTQ = gsl_matrix_alloc(m, m);
    gsl_matrix *QR = gsl_matrix_alloc(n, m);
    
    gen_rand_matrix(A);
    gsl_matrix_memcpy(Q, A);
    
    GS_decomp(Q, R);
    
    //Do multiplications for testing - QTQ and QR
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, Q, Q, 0, QTQ);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, Q, R, 0, QR);
    
    //Check that QTQ=I
    if(check_identity(QTQ, 1e-7, 1e-7) == 0){
        return_code = 1;
    }
    //Check that R is upper triangular
    if(return_code == 0){
        for(int j = 0; j<m; j++){
            for(int i = j+1; i<m; i++){
                if(equals(0, gsl_matrix_get(R, i, j), 1e-7, 1e-7) == 0){
                    return_code = 2;
                }
            }
        }   
    }
    //Check factorization, ie. that QR=A
    if(return_code == 0 && matrix_equals(A, QR, 1e-7, 1e-7) == 0){
        return_code = 3;
    }
    
    if(print == 1){
        printf("\nProduct QTQ - which should be identity\n");
        print_matrix(QTQ);
        printf("\nR matrix - which should be upper-triangular\n");
        print_matrix(R);
        printf("\nProduct QR - which of course should be equal to A\n");
        print_matrix(QR);
        printf("\nA matrix:\n");
        print_matrix(A);
    }
    
    gsl_matrix_free(A);
    gsl_matrix_free(R);
    gsl_matrix_free(Q);
    gsl_matrix_free(QR);
    gsl_matrix_free(QTQ);
    return return_code;
}
int test_partA2(int n, int print){
    int return_code = 0;
    gsl_matrix *A = gsl_matrix_alloc(n, n);
    gsl_matrix *R = gsl_matrix_calloc(n, n);
    gsl_matrix *Q = gsl_matrix_alloc(n, n);
    gsl_vector *b = gsl_vector_alloc(n);
    gsl_vector *x = gsl_vector_calloc(n);
    gsl_vector *Ax = gsl_vector_alloc(n);
    
    gen_rand_matrix(A);
    gen_rand_vector(b);
    gsl_matrix_memcpy(Q, A);
    GS_decomp(Q, R);
    
    GS_solve(Q, R, b, x);    
    gsl_blas_dgemv(CblasNoTrans, 1, A, x, 0, Ax);
    
    if(vector_equals(b, Ax, 1e-7, 1e-7) == 0){
        return_code = 1;
    }
    if(print==1){
        printf("\nChecking results of GS_solve for Ax=b - A is a new random matrix\nVector b is\n");
        print_vector(b);
        printf("\nQRx = Ax is\n");
        print_vector(Ax);
    }
    
    gsl_matrix_free(A);
    gsl_matrix_free(R);
    gsl_matrix_free(Q);
    gsl_vector_free(b);
    gsl_vector_free(x);
    gsl_vector_free(Ax);
    return return_code;
}
int test_partB(int n, int print){
    int return_code = 0;
    gsl_matrix *A = gsl_matrix_alloc(n, n);
    gsl_matrix *R = gsl_matrix_calloc(n, n);
    gsl_matrix *Q = gsl_matrix_alloc(n, n);
    gsl_matrix *B = gsl_matrix_alloc(n, n);
    gsl_matrix *AB = gsl_matrix_alloc(n, n);
    gsl_matrix *BA = gsl_matrix_alloc(n, n);
    
    gen_rand_matrix(A);
    gsl_matrix_memcpy(Q, A);
    GS_decomp(Q, R);

    GS_inverse(Q, R, B);
    
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, A, B, 0, AB);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, B, A, 0, BA);
    if(check_identity(AB, 1e-7, 1e-7) == 0){ //If AB != I return 1
        return_code = 1;
    }
    if(return_code == 0 && check_identity(BA, 1e-7, 1e-7) == 0){ //If BA != I return 2
        return_code = 2;
    }
    
    if(print == 1){
        printf("\nChecking results of GS_inverse - A is a new random matrix, B is computed inverse\nAB =\n");
        print_matrix(AB);
        printf("\nBA = \n");
        print_matrix(BA);    
    }
    
    gsl_matrix_free(A);
    gsl_matrix_free(R);
    gsl_matrix_free(Q);
    gsl_matrix_free(B);
    gsl_matrix_free(AB);
    gsl_matrix_free(BA);
    return return_code;   
}
int main(){
    srand(time(NULL));
    //A: SOLVING LINEAR EQUATIONS USING QR-DECOMPOSITION BY MOD. GS-ORTHOGONALIZATION
    printf("------------- PART A1 -------------\n\n");
    test_partA1(10, 5, 1); //n=10, m=5, print=true
    
    printf("\nRunning tests - with random m = [50, 99] and n = [m+1, m+50]:\n");
    for(int i = 0; i<10; i++){
        int m = rand() % 50 + 50;
        int n = m + rand() % 50 + 1;
        int ret = test_partA1(n, m, 0);
        printf("Test %d (size of A is %d x %d):  \t%s (code %d)\n", i+1, n, m, ret==0 ? "Success" : "Failure", ret);
    }
    printf("\n\n------------- PART A2 -------------\n\n");
    test_partA2(10, 1); // n = 10, print=true
    
    printf("\nRunning tests - with random n = [50, 200]:\n");
    for(int i = 0; i<10; i++){
        int n = rand() % 151 + 50;
        int ret = test_partA2(n, 0);
        printf("Test %d (size of A is %d x %d):  \t%s (code %d)\n", i+1, n, n, ret==0 ? "Success" : "Failure", ret);
    }
    
    //B: Matrix inverse by Gram-Schmidt QR factorization
    printf("\n\n------------- PART B  -------------\n\n");
    test_partB(10, 1);
    
    printf("\nRunning tests - with random n = [50, 200]:\n");
    for(int i = 0; i<10; i++){
        int n = rand() % 151 + 50;
        int ret = test_partB(n, 0);
        printf("Test %d (size of A is %d x %d):  \t%s (code %d)\n", i+1, n, n, ret==0 ? "Success" : "Failure", ret);
    }
    
    printf("\n\n------------- PART C  -------------\n\nTiming output sent to timing.out\n");
    FILE *toutput = fopen("timing.out", "w");
    
    for(int i = 20; i<=80; i+=2){
        int n = 10*i;
        gsl_matrix *A = gsl_matrix_alloc(n, n);
        gsl_matrix *R = gsl_matrix_alloc(n, n);
        gsl_matrix *A_copy = gsl_matrix_alloc(n,n);
        gsl_vector *tau = gsl_vector_alloc(n);
        
        gen_rand_matrix(A);
        gsl_matrix_memcpy(A_copy, A);
        
        clock_t start = clock();
        GS_decomp(A, R);
        clock_t end = clock();
        int msec_my_algo = (double)(end-start) * 1000. / CLOCKS_PER_SEC;
        
        
        start = clock();
        gsl_linalg_QR_decomp(A_copy, tau);
        end = clock();
        int msec_gsl = (double)(end-start) * 1000. / CLOCKS_PER_SEC;
        
        
        fprintf(toutput, "%i %d %d\n", n, msec_my_algo, msec_gsl);
        
        gsl_matrix_free(A);
        gsl_matrix_free(A_copy);
        gsl_matrix_free(R);
        gsl_vector_free(tau);
    }
    
    fclose(toutput);    
    return 0;
}