#include<stdio.h>
#include<math.h>
#include<assert.h>
#include<stdlib.h>
#include"spline.h"

void read_points(FILE *in_stream, int *n, double **xs, double **ys){
    int cread = fscanf(in_stream, "%d", n);
    assert(cread>0);
    assert(*n > 1);
    fprintf(stderr, "readpoints: %d\n", *n);
    
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
    int N_eval = 200;
    
    //LINEAR SPLINE
    FILE *linput = fopen("lpoints.in", "r");
    FILE *loutput = fopen("linear.out", "w");
    //Parse data from file
    int l_npoints;
    double *xs, *ys;
    read_points(linput, &l_npoints, &xs, &ys);
    
    //Evaluation of spline
    for(int i = 0; i<N_eval+1; i++){
        double z = xs[0] + (xs[l_npoints-1]-xs[0])/N_eval*i;
        double spline_y = linterp(l_npoints, xs, ys, z);
        double integrated = linterp_integ(l_npoints, xs, ys, z);
        fprintf(loutput, "%g %g %g\n", z, spline_y, integrated);
    }
    fclose(linput);
    fclose(loutput);
    free(xs);
    free(ys);
    
    //QUADRATIC SPLINE
    FILE *qinput = fopen("qpoints.in", "r");
    FILE *qoutput = fopen("quad.out", "w");
    FILE *q_coef_out = fopen("qcoef.out", "w");
    //Parse data from file
    int q_n;
    read_points(qinput, &q_n, &xs, &ys);
    qspline *s = qspline_init(q_n, xs, ys);
        
    for(int i = 0; i<N_eval+1; i++){
        double z = xs[0] + (xs[q_n-1]-xs[0])/N_eval*i;
        double spline_y = qspline_eval(s, z);
        double derivative = qspline_deriv(s, z);
        double integrated = qspline_integ(s, z);
        fprintf(qoutput, "%g %g %g %g\n", z, spline_y, derivative, integrated);
    }
    for(int i = 0; i<q_n-1; i++){ 
        fprintf(q_coef_out, "c%d = %g\n", i, s->c[i]);
        fprintf(q_coef_out, "b%d = %g\n", i, s->b[i]);  
    }
    
    qspline_free(s);
    fclose(qinput);
    fclose(qoutput);
    fclose(q_coef_out);
    free(xs);
    free(ys);
    
    return 0;
}