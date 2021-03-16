#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<stdio.h>
#include<assert.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

void print_vector(gsl_vector *v){
    int n = v->size;
    for(int i = 0; i<n; i++){
        printf("%10.7f ", gsl_vector_get(v, i));
    }
    printf("\n");
}
void print_matrix(gsl_matrix *A){
    int n = A->size1;
    int m = A->size2;
    for(int i = 0; i<n; i++){
        gsl_vector *v = gsl_vector_alloc(m);
        gsl_matrix_get_row(v, A, i);
        print_vector(v);
        gsl_vector_free(v);
    }
}
void gen_rand_vector(gsl_vector *v){
    for(int i = 0; i<v->size; i++){
        gsl_vector_set(v, i, 1.0*rand()/RAND_MAX);
    }
}
void gen_rand_matrix(gsl_matrix *A){
    int n = A->size1;
    int m = A->size2;
    for(int i = 0; i<n; i++){
        for(int j = 0; j<m; j++){
            gsl_matrix_set(A, i, j, 1.0*rand()/RAND_MAX);
        }
    }
}
void GS_decomp(gsl_matrix *A, gsl_matrix *R){
    int n = A->size1;
    int m = A->size2;
    assert(R->size1==R->size2 && R->size1 == m);
    assert(n>=m);
    for(int i = 0; i<m; i++){
        gsl_vector *a_i = gsl_vector_alloc(n);
        gsl_vector *q_i = gsl_vector_calloc(n);
        
        gsl_matrix_get_col(a_i, A, i);
        
        double a_i_norm = gsl_blas_dnrm2(a_i);
        gsl_matrix_set(R, i, i, a_i_norm);
        
        gsl_blas_daxpy(1./a_i_norm, a_i, q_i);
        gsl_matrix_set_col(A, i, q_i);
        
        for(int j = i+1; j<m; j++){
            gsl_vector *a_j = gsl_vector_alloc(n);
            gsl_matrix_get_col(a_j, A, j);
            
            double inner_prod_qa = 0;
            gsl_blas_ddot(q_i, a_j, &inner_prod_qa);
            gsl_matrix_set(R, i, j, inner_prod_qa);
            
            gsl_blas_daxpy(-inner_prod_qa, q_i, a_j);
            gsl_matrix_set_col(A, j, a_j);
            gsl_vector_free(a_j);
        }
        gsl_vector_free(a_i);
        gsl_vector_free(q_i);
    }
}
void backsub(gsl_matrix *U, gsl_vector *c){
    assert(U->size1 == U->size2);
    int n = U->size1;
    for(int i = n-1; i>=0; i--){
        double s = gsl_vector_get(c, i);
        for(int k = i+1; k<n; k++){
            s -= gsl_matrix_get(U, i, k)*gsl_vector_get(c, k);
        }
        gsl_vector_set(c, i, s/gsl_matrix_get(U,i,i));
    }
}
void GS_solve(gsl_matrix *Q, gsl_matrix *R, gsl_vector *b, gsl_vector *x){
    gsl_blas_dgemv(CblasTrans, 1, Q, b, 0, x);
    backsub(R, x);
}
int equals(double a, double b, double tau, double epsilon){
    if(fabs(a-b)<tau){
        return 1;
    }
    if(fabs(a-b)/(fabs(a)+fabs(b)) < epsilon/2){
        return 1;
    }
    return 0;
}
int matrix_equals(gsl_matrix *A, gsl_matrix *B, double tau, double eps){
    int n1 = A->size1; int m1 = A->size2;
    int n2 = B->size1; int m2 = B->size2;
    assert(n1 == n2 && m1 == m2);
    for(int i = 0; i<n1; i++){
        for(int j = 0; j<m1; j++){
            double aij = gsl_matrix_get(A, i, j);
            double bij = gsl_matrix_get(B, i, j);
            if(equals(aij, bij, tau, eps) == 0){
                return 0;
            }
        }
    }
    return 1;
}
int vector_equals(gsl_vector *v, gsl_vector *u, double tau, double eps){
    int n1 = v->size; int n2 = u->size;
    assert(n1 == n2);
    for(int i = 0; i<n1; i++){
        double vi = gsl_vector_get(v, i);
        double ui = gsl_vector_get(u, i);
        if(equals(vi, ui, tau, eps) == 0){
            return 0;
        }
    }
    return 1;
}
void init_identity(gsl_matrix *I){
    assert(I->size1 == I->size2);
    int n = I->size1;
    for(int i = 0; i<n; i++){
        gsl_matrix_set(I, i, i, 1);
    }
}
//Returns 0 on success
//Rerturns non-zero and different for each possible failure
int test_partA1(int seed){
    srand(seed);
    int m=rand() % 51 + 10;
    int n = m + rand() % 50 + 1;
    gsl_matrix *A = gsl_matrix_alloc(n, m);
    gsl_matrix *R = gsl_matrix_calloc(m, m);
    gsl_matrix *Q = gsl_matrix_alloc(n, m);
    gsl_matrix *I = gsl_matrix_calloc(m, m);
    gsl_matrix *QTQ = gsl_matrix_alloc(m, m);
    gsl_matrix *QR = gsl_matrix_alloc(n, m);
    
    init_identity(I);
    gen_rand_matrix(A);
    gsl_matrix_memcpy(Q, A);
    
    GS_decomp(Q, R);
    
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, Q, Q, 0, QTQ);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, Q, R, 0, QR);
    if(matrix_equals(I, QTQ, 1e-7, 1e-7) == 0){
        return 1;
    }
    
    //TODO: Check than R is upper triangular
    
    if(matrix_equals(A, QR, 1e-7, 1e-7) == 0){
        return 2;
    }
    
    gsl_matrix_free(A);
    gsl_matrix_free(R);
    gsl_matrix_free(Q);
    gsl_matrix_free(I);
    gsl_matrix_free(QR);
    gsl_matrix_free(QTQ);
    return 0;
}
int test_partA2(int seed){
    srand(seed);
    int n = rand() % 90 + 10;
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
        return 1;
    }
    
    gsl_matrix_free(A);
    gsl_matrix_free(R);
    gsl_matrix_free(Q);
    gsl_vector_free(b);
    gsl_vector_free(x);
    gsl_vector_free(Ax);
    return 0;
}
int main(){
    srand(time(NULL));
    //A: SOLVING LINEAR EQUATIONS USING QR-DECOMPOSITION BY MOD. GS-ORTHOGONALIZATION
    {
        printf("------------- PART A1 -------------\n\n");
        int n=10,m=5;
        gsl_matrix *A = gsl_matrix_alloc(n, m);
        gsl_matrix *R = gsl_matrix_calloc(m, m);
        gsl_matrix *Q = gsl_matrix_alloc(n, m);
        gsl_matrix *QTQ = gsl_matrix_alloc(m, m);
        gsl_matrix *QR = gsl_matrix_alloc(n, m);
        
        gen_rand_matrix(A);
        gsl_matrix_memcpy(Q, A);
        
        GS_decomp(Q, R);
        
        gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, Q, Q, 0, QTQ); //Calculate QTQ (=I)
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, Q, R, 0, QR); //Calculate QR (=A)
        
        //Printing results
        printf("\nProduct QTQ - which should be identity\n");
        print_matrix(QTQ);
        printf("\nR matrix - which should be upper-triangular\n");
        print_matrix(R);
        printf("\nProduct QR - which of course should be equal to A\n");
        print_matrix(QR);
        printf("\nA matrix:\n");
        print_matrix(A);
        
        gsl_matrix_free(A);
        gsl_matrix_free(Q);
        gsl_matrix_free(R);
        gsl_matrix_free(QR);
        gsl_matrix_free(QTQ);
        
        printf("\nRunning tests - with random m = [10, 60] and n = [m+1, m+50]:\n");
        for(int i = 0; i<10; i++){
            int ret = test_partA1(i+time(NULL));
            printf("Test %d: %s (code %d)\n", i+1, ret==0 ? "Success" : "Failure", ret);
        }
    }
    {
        printf("\n\n------------- PART A2 -------------\n\n");
        
        int n = 10;
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
        gsl_blas_dgemv(CblasNoTrans, 1, A, x, 0, Ax); //Check by calculating Ax (=b)
        
        printf("\nChecking results of GS_solve for Ax=b - A is a new random matrix\nVector b is\n");
        print_vector(b);
        printf("\nQRx = Ax is\n");
        print_vector(Ax);
        
        gsl_vector_free(b);
        gsl_vector_free(x);
        gsl_vector_free(Ax);
        gsl_matrix_free(A);
        gsl_matrix_free(R);
        gsl_matrix_free(Q);
        
        printf("\nRunning tests - with random n = [10, 99]:\n");
        for(int i = 0; i<10; i++){
            int ret = test_partA2(i + time(NULL)+42);
            printf("Test %d: %s (code %d)\n", i+1, ret==0 ? "Success" : "Failure", ret);
        }
    }
    return 0;
}