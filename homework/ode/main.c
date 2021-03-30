#include<stdio.h>
#include"ode.h"
#include"matrix.h"

void func_shm(double t, gsl_vector *y, gsl_vector *dydt){
    //u''=-u
    //y0 = u
    //y1 = u'
    //
    //y0' = y1
    //y1' = -y0
    
    gsl_vector_set(dydt, 0, gsl_vector_get(y, 1));
    gsl_vector_set(dydt, 1, -gsl_vector_get(y, 0));
}

int main(){
    int n = 4000, m=2;
    // int m = 2;
    gsl_matrix *ylist = gsl_matrix_alloc(n, m);
    gsl_vector *xlist = gsl_vector_alloc(n);
    gsl_vector *yt = gsl_vector_alloc(m);
    gsl_vector *yh = gsl_vector_alloc(m);
    gsl_vector *err = gsl_vector_alloc(m);
    
    //Initial value
    gsl_matrix_set(ylist, 0, 0, 1);
    gsl_matrix_set(ylist, 0, 1, 0);
    // gsl_vector_set(yt, 0, 1);
    // gsl_vector_set(yt, 1, 0);
    
    // FILE *data_out = fopen("data.out", "w");
    // double t = 0;
    // for(int i = 0; i<100; i++){
    //     rkstep45(func_shm, t, yt, 0.002, yh, err);
    //     fprintf(data_out, "%g %g %g\n", t, gsl_vector_get(yt, 0), gsl_vector_get(yt, 1));
    //     gsl_vector_memcpy(yt, yh);
    //     t += 0.01;
    // }
    int k = driver(func_shm, 0, 20, 0.001, 1e-3, 1e-5, ylist, xlist);
    FILE *data_out = fopen("data.out", "w");
    for(int i = 0; i<abs(k); i++){
        gsl_vector_view yi = gsl_matrix_row(ylist, i);
        fprintf(data_out, "%g %g %g\n", gsl_vector_get(xlist, i), gsl_vector_get(&yi.vector, 0), gsl_vector_get(&yi.vector, 1));
    }
    gsl_vector_free(yt);
    gsl_vector_free(yh);
    gsl_vector_free(err);
    fclose(data_out);
    return 0;
}