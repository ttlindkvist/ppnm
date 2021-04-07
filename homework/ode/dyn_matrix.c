#include "dyn_matrix.h"
#include<stdlib.h>
#include<stdio.h>
#include<assert.h>

//Dynamic matrix implementation
dyn_matrix *dyn_matrix_alloc(unsigned int n1, unsigned int n2){
    dyn_matrix *A = (dyn_matrix*)malloc(sizeof(dyn_matrix));
    A->n1 = n1;
    A->n2 = n2;
    A->data = (double*)malloc(sizeof(double)*n1*n2);
    return A;
}
void dyn_matrix_add_rows(dyn_matrix *A, unsigned int n){
    A->data = (double*)realloc(A->data, sizeof(double)*(A->n1+n)*A->n2);
    A->n1 += n;
}
void dyn_matrix_set(dyn_matrix *A, unsigned int i, unsigned int j, double x){
    assert(i < A->n1 && j < A->n2);
    A->data[i*A->n2+j] = x;
}
double dyn_matrix_get(dyn_matrix *A, unsigned int i, unsigned int j){
    assert(i < A->n1 && j < A->n2);
    return A->data[i*A->n2+j];
}

double *dyn_matrix_row(dyn_matrix *A, unsigned int i){
    return A->data + i*A->n2;
}

void dyn_matrix_free(dyn_matrix *A){
    free(A->data);
    free(A);
}
void dyn_matrix_print(dyn_matrix *A){
    for(int i = 0; i<A->n1; i++){
        for(int j = 0; j<A->n2; j++){
            printf("%g ", A->data[i*A->n2+j]);
        }
        printf("\n");
    }
}