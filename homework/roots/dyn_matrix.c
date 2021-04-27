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
gsl_vector_view dyn_matrix_row_view(dyn_matrix *A, unsigned int i){
    assert(i<A->n1);
    gsl_vector_view view = gsl_vector_view_array(A->data + i*A->n2, A->n2);
    return view;
}

void dyn_matrix_free(dyn_matrix *A){
    free(A->data);
    free(A);
}
void dyn_matrix_print(dyn_matrix *A){
    for(unsigned int i = 0; i<A->n1; i++){
        for(unsigned int j = 0; j<A->n2; j++){
            printf("%g ", A->data[i*A->n2+j]);
        }
        printf("\n");
    }
}

dyn_vector *dyn_vector_alloc(unsigned int n){
    dyn_vector *v = (dyn_vector*)malloc(sizeof(dyn_vector));
    v->n = n;
    v->data = (double*)malloc(sizeof(double)*n);
    return v;
}
void dyn_vector_free(dyn_vector *v){
    free(v->data);
    free(v);
}

void dyn_vector_inc_size(dyn_vector *v, unsigned int n){
        v->data = (double*)realloc(v->data, sizeof(double)*(v->n+n));
        v->n += n;
}
double dyn_vector_get(dyn_vector *v, unsigned int i){
    assert(i<v->n);
    return v->data[i];
}
void dyn_vector_set(dyn_vector *v, unsigned int i, double x){
    assert(i<v->n);
    v->data[i] = x;
}
gsl_vector_view dyn_vector_to_vec_view(dyn_vector *v){
    gsl_vector_view view = gsl_vector_view_array(v->data, v->n);
    return view;
}