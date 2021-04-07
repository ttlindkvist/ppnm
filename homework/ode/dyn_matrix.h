#ifndef __DYN_MATRIX_H__
#define __DYN_MATRIX_H__

typedef struct {int n1; int n2; double* data;} dyn_matrix;

dyn_matrix *dyn_matrix_alloc(unsigned int n1, unsigned int n2);
void dyn_matrix_add_rows(dyn_matrix *A, unsigned int n);
double dyn_matrix_get(dyn_matrix *A, unsigned int i, unsigned int j);
void dyn_matrix_set(dyn_matrix *A, unsigned int i, unsigned int j, double x);

double *dyn_matrix_row(dyn_matrix *A, unsigned int i);

void dyn_matrix_free(dyn_matrix *A);

void dyn_matrix_print(dyn_matrix *A);

#endif