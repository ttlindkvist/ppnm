#include"matrix.h"
#include<gsl/gsl_blas.h>
#include<gsl/gsl_eigen.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<time.h>

//Checks that the matrix A is diagnonal.
//It is known that A is symmetric, so only upper triangle is checked
int checkDiagonal(gsl_matrix *A){
    assert(A->size1 == A->size2);
    int n = A->size1;
    //We know that A is symmetric - only check upper triangle
    for(int i = 0; i<n; i++){
        for(int j = i+1; j<n; j++){
            // Quite large epsilon, but we know that further rotation doesn't change the diagonal/calculated eigenvalues
            if(equals(gsl_matrix_get(A, i, j), 0, 1e-6, 1e-6) == 0){
                printf("Non-zeroed element (%d,%d) is in fact: %g\n", i, j, gsl_matrix_get(A, i, j));
                return 0;
            }
        }
    }
    return 1;
}

int testA(int n, int print){
    gsl_matrix *A = gsl_matrix_alloc(n, n);
    gsl_matrix *D = gsl_matrix_alloc(n, n);
    gsl_matrix *V = gsl_matrix_alloc(n, n);
    gsl_matrix *temp = gsl_matrix_alloc(n, n);
    gsl_matrix *VTV = gsl_matrix_alloc(n, n);
    gsl_matrix *VTAV = gsl_matrix_alloc(n, n);
    gsl_matrix *VDVT = gsl_matrix_alloc(n, n);
    gen_rand_symm_matrix(A);
    gsl_matrix_memcpy(D, A);
    
    jacobi_diag(D, V);
    if(print==1) print_matrix(D);
    
    //Check that D is diagonal
    if(checkDiagonal(D)==0){
        return 1;
    }
    //Check that VTAV=D
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, V, A, 0, temp);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, temp, V, 0, VTAV);
    if(matrix_equals(VTAV, D, 1e-7, 1e-7)==0){
        return 2;
    }
    //Check that VDVT=A
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, V, D, 0, temp);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, temp, V, 0, VDVT);
    if(matrix_equals(VDVT, A, 1e-7, 1e-7)==0){
        return 3;
    }
    //Check that VTV=I
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, V, V, 0, VTV);
    if(check_identity(VTV, 1e-7, 1e-7)==0){
        return 4;
    }
    
    if(print == 1){
        printf("\nPerforming eigenvalue decomposition on A\nA is\n");
        print_matrix(A);
        printf("\nProduct VDVT is A\n");
        print_matrix(VDVT);
        printf("\nD is diagonal\n");
        print_matrix(D);
        printf("\nProduct VTAV is D\n");
        print_matrix(VTAV);
        printf("\nAnd VTV is identity\n");
        print_matrix(VTV);
    }
    
    gsl_matrix_free(A);
    gsl_matrix_free(D);
    gsl_matrix_free(V);
    gsl_matrix_free(temp);
    gsl_matrix_free(VTV);
    gsl_matrix_free(VDVT);
    gsl_matrix_free(VTAV);
    return 0;
}
void partB(FILE *eigenfunc_out){
    int n = 200;
    double s=1.0/(n+1);
    gsl_matrix *H = gsl_matrix_calloc(n, n);
    gsl_matrix *V = gsl_matrix_alloc(n, n);
    
    //Initialize hamiltonian matrix
    for(int i = 0; i<n-1; i++){
        gsl_matrix_set(H,i,i,-2);
        gsl_matrix_set(H,i,i+1,1);
        gsl_matrix_set(H,i+1,i,1);
    }
    gsl_matrix_set(H,n-1,n-1,-2);
    gsl_matrix_scale(H,-1/s/s);
    
    jacobi_diag(H, V);
        
    printf("Diagonalized hamiltonian of size %d\nThe lowest 5 eigenenergies are\n", n);
    for (int k = 0; k < 5; k++){
        double exact = M_PI*M_PI*(k+1)*(k+1); //pi²*n² in these units
        double calculated = gsl_matrix_get(H,k,k);
        printf("n=%i \tcalculated to be = %.10g \texact = %.10g\n", k+1, calculated, exact);
    }
    
    int signs[3] = {-1, -1, 1}; //some signs of eigenfunctions needs to be flipped
    for(int k = 0; k<3; k++){
        fprintf(eigenfunc_out, "#INDEX %i - computed eigenfunction %i\n", k, k);
        fprintf(eigenfunc_out, "0 0\n");
        for(int i=0;i<n;i++){
            double x = (i+1.0)/(n+1);
            double calc = signs[k]*1/sqrt(s)*gsl_matrix_get(V,i,k); // 1/sqrt(s) is normalization factor
            fprintf(eigenfunc_out, "%g %g\n", x, calc);            
        }
        fprintf(eigenfunc_out, "1 0\n\n\n");
    }
    
    for(int k = 0; k<3; k++){
        fprintf(eigenfunc_out, "#INDEX %i - analytic eigenfunction %i\n", k+3, k);
        for(int i=0;i<=1000;i++){
            double x = i/1000.;
            double exact = sqrt(2)*sin((k+1)*M_PI*x);
            fprintf(eigenfunc_out, "%g %g\n", x, exact);            
        }
        fprintf(eigenfunc_out, "\n\n\n");
    }
    printf("\nPlot produced - eigenfunc.out\n\n");
    gsl_matrix_free(H);
    gsl_matrix_free(V);
}
int main(){
    srand(time(NULL));
    printf("------------- PART A -------------\n\n");
    printf("Using non-chaninging diagonal as convergence-criteria\n\n");
    fprintf(stderr, "Running part A\n");
    testA(5, 1);

    printf("\nPerforming tests on part A - with n=[100, 300]\n");
    for(int i = 0; i<5; i++){
        int n = rand() % 201 + 100;
        fprintf(stderr, "Running test %d with (%d x %d)\n", i+1, n, n);
        int ret = testA(n, 0);
        printf("Test %d (size of A is %d x %d):  \t%s (code %d)\n", i+1, n, n, ret==0 ? "Success" : "Failure", ret);
    }
    
    printf("\n\n------------- PART B -------------\n\n");
    FILE *eigenfunc_out = fopen("eigenfunc.out", "w");
    fprintf(stderr, "Running part B\n");
    partB(eigenfunc_out);
    fclose(eigenfunc_out);
    
    printf("\n\n------------- PART C -------------\n\n");
    fprintf(stderr, "Running part C\n");
    //Time my implementation vs gsl's
    FILE *timing_out = fopen("timing.out", "w");
    for(int i = 1; i<=20; i++){
        int n = i*10;
        
        
        gsl_matrix *A = gsl_matrix_alloc(n, n);
        gsl_matrix *A_copy = gsl_matrix_alloc(n, n);
        gsl_matrix *A_copy2 = gsl_matrix_alloc(n, n);
        gsl_matrix *A_copy3 = gsl_matrix_alloc(n, n);
        gsl_matrix *V = gsl_matrix_alloc(n, n);
        
        gen_rand_symm_matrix(A);
        gsl_matrix_memcpy(A_copy, A);
        gsl_matrix_memcpy(A_copy2, A);
        gsl_matrix_memcpy(A_copy3, A);
        
        clock_t t = clock();
        jacobi_diag(A, V);
        double j_duration_ms = (double)(clock()-t)/1000.;
        
        t = clock();
        jacobi_diag_sum(A_copy2, V, 1e-6);
        double j_sum_duration_ms = (double)(clock()-t)/1000.;
        
        t = clock();
        jacobi_diag_opt(A_copy3, V);
        double j_opt_duration_ms = (double)(clock()-t)/1000.;
        
        gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(n); 
        gsl_vector *eval = gsl_vector_alloc(n);
        
        t = clock();
        gsl_eigen_symmv(A_copy, eval, V, w);
        double gsl_duration_ms = (double)(clock() - t)/1000.;
        
        fprintf(timing_out, "%d %g %g %g %g\n", n, j_duration_ms, gsl_duration_ms, j_sum_duration_ms, j_opt_duration_ms);
        if(i%4==0){
            fprintf(stderr, "Part C progress (%d/20)\n", i);            
        }
        gsl_matrix_free(A);
        gsl_matrix_free(A_copy);
        gsl_matrix_free(A_copy2);
        gsl_matrix_free(A_copy3);
        gsl_matrix_free(V);
        
        gsl_eigen_symmv_free(w);
        gsl_vector_free(eval);
    }
    fclose(timing_out);
    printf("Timing output to file timing.out - plots generated.\n");
    printf("Plot of our Jacobi diagonalization in jacobi.timing.png, and gsl implementation in gsl.timing.png\n");
    printf("As seen from the plots the gsl implementation is much faster\n");
    //Optimize my implementation
    return 0;
}