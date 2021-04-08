#ifndef __DYN_MATRIX_H__
#define __DYN_MATRIX_H__
#include<gsl/gsl_vector.h>

typedef struct {int n1; int n2; double* data;} dyn_matrix;
typedef struct {int n; double* data;} dyn_vector;

dyn_matrix *dyn_matrix_alloc(unsigned int n1, unsigned int n2);
void dyn_matrix_add_rows(dyn_matrix *A, unsigned int n);
double dyn_matrix_get(dyn_matrix *A, unsigned int i, unsigned int j);
void dyn_matrix_set(dyn_matrix *A, unsigned int i, unsigned int j, double x);

double *dyn_matrix_row(dyn_matrix *A, unsigned int i);
gsl_vector_view dyn_matrix_row_view(dyn_matrix *A, unsigned int i);

void dyn_matrix_free(dyn_matrix *A);

void dyn_matrix_print(dyn_matrix *A);

dyn_vector *dyn_vector_alloc(unsigned int n);
void dyn_vector_free(dyn_vector *v);
void dyn_vector_inc_size(dyn_vector *v, unsigned int n);
double dyn_vector_get(dyn_vector *v, unsigned int i);
void dyn_vector_set(dyn_vector *v, unsigned int i, double x);
gsl_vector_view dyn_vector_to_vec_view(dyn_vector *v);

#endif