#include<stdio.h>
#include"dyn_matrix.h"

int main(){
    int n1=10,n2=10;
    dyn_matrix *A = dyn_matrix_alloc(n1, n2);
    for(int i = 0; i<n1; i++){
        for(int j = 0; j<n2; j++){
            dyn_matrix_set(A, i, j, i+j);
        }
    }
    dyn_matrix_add_rows(A, 1);
    for(int i = 0; i<n2; i++){
        dyn_matrix_set(A, n1, i, i+100);
    }
    dyn_matrix_print(A);
    
    printf("\nTake row 3:\n");
    double *row = dyn_matrix_row(A, 2);
    for(int i = 0; i<n2; i++){
        printf("%g ", row[i]);
    }
    printf("\n\n");
    
    
    dyn_matrix_free(A);
}