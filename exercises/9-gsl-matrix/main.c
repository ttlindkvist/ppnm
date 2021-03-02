#include<stdio.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_eigen.h>

// A: Solve the following
// [	 6.13 	 -2.90 	 5.86 	] [x0] 	 	 [6.23]
// [	 8.08 	 -6.31 	 -3.89 	] [x1] 	 = 	 [5.37]
// [	 -4.36 	 1.00 	 0.19 	] [x2] 	 	 [2.29]

void print_vector(gsl_vector *x, int n){
    for(int i = 0; i<n; i++){
        printf("%g ", gsl_vector_get(x, i));
    }
}
void init_Hilbert_matrix(gsl_matrix *a, int n){
    for(int i = 0; i<n; i++){
        for(int j = 0; j<n; j++){
            gsl_matrix_set(a, i, j, 1.0/(1+i+j));
        }
    }
}

int main(){
    double a_data[] = {
        6.13, -2.90, 5.86,
        8.08, -6.31, -3.89,
        -4.36, 1.00, 0.19
    };
    double b_data[] = {6.23, 5.37, 2.29};
    gsl_matrix_view A = gsl_matrix_view_array(a_data, 3, 3);
    gsl_vector_view b = gsl_vector_view_array(b_data, 3);
    gsl_vector *x = gsl_vector_alloc(3);
    gsl_matrix *A_copy = gsl_matrix_alloc(3, 3);
    gsl_matrix_memcpy(A_copy, &A.matrix); 
    
    gsl_linalg_HH_solve(&A.matrix, &b.vector, x);
    
    gsl_vector *product = gsl_vector_alloc(3);
    gsl_blas_dgemv(CblasNoTrans, 1, A_copy, x, 0, product);
    
    printf("x = ");
    print_vector(x, 3);
    printf("\n");
    
    printf("A*x = ");
    print_vector(product, 3);
    printf("\n");
    
    printf("b = ");
    print_vector(&b.vector, 3);
    printf("\n");
    
    
    gsl_vector_free(x);
    gsl_vector_free(product);
    gsl_matrix_free(A_copy);
    
    printf("\n4th order Hilbert matrix eigenvalues and eigenvectors\n\n");
    
    //Compute eigenvectors and eigenvalues for 4th order Hilbert matrix.
    gsl_matrix *H = gsl_matrix_alloc(4, 4);
    init_Hilbert_matrix(H, 4);
    
    gsl_vector *eval = gsl_vector_alloc(4);
    gsl_matrix *evec = gsl_matrix_alloc(4, 4);

    gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(4);

    gsl_eigen_symmv(H, eval, evec, w);

    gsl_eigen_symmv_free(w);

    gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC);

    for (int i = 0; i < 4; i++){
        double eval_i = gsl_vector_get(eval, i);
        gsl_vector_view evec_i = gsl_matrix_column(evec, i);

        printf("eigenvalue = %g\n", eval_i);
        printf("eigenvector = ");
        print_vector(&evec_i.vector, 4);
        printf("\n\n");
    }

    gsl_vector_free(eval);
    gsl_matrix_free(evec);
    gsl_matrix_free(H);
    
    return 0;
}