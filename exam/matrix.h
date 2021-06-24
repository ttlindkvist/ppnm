#ifndef __MATRIX__H__
#define __MATRIX__H__
#include<gsl/gsl_matrix.h>

void print_vector(gsl_vector *v);
void print_matrix(gsl_matrix *A);
void gen_rand_vector(gsl_vector *v);
void gen_rand_matrix(gsl_matrix *A);
int equals(double a, double b, double tau, double epsilon);
int vector_equals(gsl_vector *v, gsl_vector *u, double tau, double eps);
int matrix_equals(gsl_matrix *A, gsl_matrix *B, double tau, double eps);
int check_identity(gsl_matrix *A, double tau, double eps);
void backsub(gsl_matrix *U, gsl_vector *c);
void GS_decomp(gsl_matrix *A, gsl_matrix *R);
void GS_solve(gsl_matrix *Q, gsl_matrix *R, gsl_vector *b, gsl_vector *x);
void GS_inverse(gsl_matrix *Q, gsl_matrix *R, gsl_matrix *B);

//Jacobi diagonalization
void gen_rand_symm_matrix(gsl_matrix *A);
void timesJ(gsl_matrix *A, int p, int q, double theta);
void Jtimes(gsl_matrix *A, int p, int q, double theta);

int SVD_two_jaco_square(gsl_matrix *A, gsl_matrix *V, gsl_matrix *U);

#endif