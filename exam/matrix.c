#include"matrix.h"
#include<assert.h>
#include<gsl/gsl_blas.h>
#include<math.h>


void print_vector(gsl_vector *v){
    int n = v->size;
    for(int i = 0; i<n; i++){
        printf("%10.7f ", gsl_vector_get(v, i));
    }
    printf("\n");
}
void print_matrix(gsl_matrix *A){
    int n = A->size1;
    for(int i = 0; i<n; i++){
        gsl_vector_view v  = gsl_matrix_row(A, i);
        print_vector(&v.vector);
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

int equals(double a, double b, double tau, double epsilon){
    if(fabs(a-b)<tau){
        return 1;
    }
    if(fabs(a-b)/(fabs(a)+fabs(b)) < epsilon/2){
        return 1;
    }
    return 0;
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
int check_identity(gsl_matrix *A, double tau, double eps){
    assert(A->size1==A->size2);
    int n = A->size1;
    for(int i = 0; i<n; i++){
        for(int j = 0; j<n; j++){
            double aij = gsl_matrix_get(A, i, j);
            if(i==j){
                if(equals(1, aij, tau, eps) == 0){
                    return 0;
                }
            } else {
                if(equals(0, aij, tau, eps) == 0){
                    return 0;
                } 
            }
        }
    }
    return 1;
}
int check_diagonal(gsl_matrix *A, double tau, double eps){
    assert(A->size1==A->size2);
    int n = A->size1;
    for(int i = 0; i<n; i++){
        for(int j = 0; j<n; j++){
            if(i!=j){
                double aij = gsl_matrix_get(A, i, j);
                if(equals(0, aij, tau, eps) == 0){
                    return 0;
                } 
            }
        }
    }
    return 1;
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
void GS_decomp(gsl_matrix *A, gsl_matrix *R){
    int n = A->size1;
    int m = A->size2;
    assert(R->size1==R->size2 && R->size1 == m);
    assert(n>=m);
    
    //Reset R
    gsl_matrix_set_all(R, 0);
    
    for(int i = 0; i<m; i++){
        gsl_vector_view a_i = gsl_matrix_column(A, i);
        
        double a_i_norm = gsl_blas_dnrm2(&a_i.vector);
        gsl_matrix_set(R, i, i, a_i_norm);
        
        //Normalize column i in A 
        gsl_vector_scale(&a_i.vector, 1./a_i_norm);
                
        for(int j = i+1; j<m; j++){
            gsl_vector_view a_j = gsl_matrix_column(A, j);
            
            double inner_prod_qa = 0;
            gsl_blas_ddot(&a_i.vector, &a_j.vector, &inner_prod_qa);
            gsl_matrix_set(R, i, j, inner_prod_qa);
            
            // a_j <- a_j - <q_i|a_j>q_i
            gsl_blas_daxpy(-inner_prod_qa, &a_i.vector, &a_j.vector);
        }
    }
}
void GS_solve(gsl_matrix *Q, gsl_matrix *R, gsl_vector *b, gsl_vector *x){
    //Ax=(QR)x=b, QTQ=I
    //Rx=QTb
    //R is upper triangular -> use backsub
    gsl_blas_dgemv(CblasTrans, 1, Q, b, 0, x);
    backsub(R, x);
}
void GS_inverse(gsl_matrix *Q, gsl_matrix *R, gsl_matrix *B){
    assert(Q->size1 == Q->size2 && R->size1 == R->size2 && Q->size1 == R->size1);
    int n = Q->size1;
    assert(B->size1 == n && B->size2 == n);
    gsl_vector *e_i = gsl_vector_alloc(n);
    gsl_vector *x = gsl_vector_alloc(n);
    for(int i = 0; i<n; i++){ // Solve Ax_i = e_i where A⁻¹ = {x_i}, A = QR
        gsl_vector_set_basis(e_i, i);
        gsl_vector_set_zero(x);
        
        gsl_vector_set(e_i, i, 1);
        GS_solve(Q, R, e_i, x);
        gsl_matrix_set_col(B, i, x);
    }
    gsl_vector_free(e_i);
    gsl_vector_free(x);
}