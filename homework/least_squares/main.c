#include"matrix.h"
#include<stdio.h>
#include<math.h>
#include<gsl/gsl_blas.h>
#include<assert.h>

double f0(double x){
    return 1.;
}
double f1(double x){
    return x;
}

void LS_fit(double *x, double *y, double *yerr, int n, double (**fs)(double), int m, gsl_vector *c, gsl_matrix *S){
    // printf("%g\n", fs[0](1));
    gsl_matrix *A = gsl_matrix_alloc(n, m);
    gsl_matrix *R = gsl_matrix_alloc(m, m);
    gsl_vector *b = gsl_vector_alloc(n);
    
    for(int i = 0; i<n; i++){
        for(int k = 0; k<m; k++){
            gsl_matrix_set(A, i, k, fs[k](x[i])/yerr[i]);
        }
        gsl_vector_set(b, i, y[i]/yerr[i]);
    }
    //Decompose A into Q(in A's place) and R 
    GS_decomp(A, R);
    //Solve Rc=QTb for c
    GS_solve(A, R, b, c);
    
    //Find covariance matrix S = R⁻¹(R⁻¹)^T
    //First R inverse
    gsl_matrix *R_inv = gsl_matrix_alloc(m, m);
    gsl_matrix *I = gsl_matrix_alloc(m, m);
    gsl_matrix_set_identity(I);
    GS_inverse(I, R, R_inv);
    
    //Now compute product
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, R_inv, R_inv, 0, S);
    
    gsl_matrix_free(A);
    gsl_matrix_free(R);
    gsl_matrix_free(I);
    gsl_vector_free(b);
}

void read_points(FILE *in_stream, int *n, double **xs, double **ys){
    int cread = fscanf(in_stream, "%d", n);
    assert(cread>0);
    assert(*n > 1);
    
    *xs = (double*)malloc((*n)*sizeof(double));
    *ys = (double*)malloc((*n)*sizeof(double));
    
    double x, y;
    int i = 0;
    while(fscanf(in_stream, "%lg %lg", &x, &y) > 0){
        assert(i<*n);        
        (*xs)[i] = x;
        (*ys)[i] = y;
        i++;
    }
    assert(i==*n);
}

int main(){
    FILE *input = fopen("data.in", "r");
    
    int n;
    double *xs;
    double *ys;
    printf("\nLoading datapoints from data.in\n");
    read_points(input, &n, &xs, &ys);
    
    fclose(input);

    double ys_log[n];
    double yerr_log[n];
    for(int i = 0; i<n; i++){
        ys_log[i] = log(ys[i]);
        yerr_log[i] = 1./20; //(y[i]/20)/y[i]
    }
    
    const int m = 2;
    double (*fs[])(double) = {f0, f1};
    gsl_vector *c = gsl_vector_alloc(m);
    gsl_matrix *S = gsl_matrix_alloc(m, m);
    
    LS_fit(xs, ys_log, yerr_log, n, fs, m, c, S);
    
    printf("\nFit data to ln(y) = c_0 + c_1*x => y = exp(c_0)*exp(c_1*x) = a*exp(-gamma*x)\n");
    printf("Coefficient vector c\n");
    print_vector(c);
    
    double a = exp(gsl_vector_get(c, 0));
    double a_err = a*sqrt(gsl_matrix_get(S, 0, 0));
    double g = -gsl_vector_get(c, 1);
    double g_err = sqrt(gsl_matrix_get(S, 1, 1));
    
    printf("\nThis gives the function: f(x) = %g*exp(-%g*x)\n", a, g);
    
    printf("\nCovariance matrix S\n");
    print_matrix(S);
    
    printf("\nThis gives a relative activity @ time 0: a = %g +- %g\n", a, a_err);
    printf("And gamma = %g +- %g\n", g, g_err);
    printf("Which gives a half-life of (%g +- %g) d\n", log(2)/g, log(2)/(g*g)*g_err);
    printf("Half-life for 224Ra = 3.63 d\nWhich is not terribly far away from our fitted value.\nHowever it still outside the estimated uncertainty.\n");
    gsl_vector_free(c);
    gsl_matrix_free(S);
    free(xs);
    free(ys);
    
    //Make output file for plotting
    FILE *output = fopen("fit.out", "w");
    int N = 1000;
    fprintf(output, "#index 0 - best fit\n");
    for(int i = 0; i<N; i++){
        double x = 16./N * i;
        fprintf(output, "%g %g\n", x, a*exp(-g*x));
    }
    fprintf(output, "\n\n#index 1 - fit + deltac\n");
    for(int i = 0; i<N; i++){
        double x = 16./N * i;
        fprintf(output, "%g %g\n", x, (a+a_err)*exp(-(g-g_err)*x));
    }
    fprintf(output, "\n\n#index 2 - fit - deltac\n");
    for(int i = 0; i<N; i++){
        double x = 16./N * i;
        fprintf(output, "%g %g\n", x, (a-a_err)*exp(-(g+g_err)*x));
    }
    fclose(output);
    printf("From fig_plot.png the points is seen together with the range specified by (c+-delta c)\n");
    return 0;
}