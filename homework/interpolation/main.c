#include<stdio.h>
#include<math.h>
#include<assert.h>
#include<stdlib.h>
#include"spline.h"
#include<gsl/gsl_interp.h>

//Read points from file open in FILE* in_stream
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
    
    //A: LINEAR SPLINE
    FILE *linput = fopen("lpoints.in", "r");
    FILE *loutput = fopen("linear.out", "w");
    //Parse data from file
    int l_n;
    double *xs, *ys;
    read_points(linput, &l_n, &xs, &ys);
    
    gsl_interp *linear = gsl_interp_alloc(gsl_interp_linear, l_n);
    gsl_interp_init(linear, xs, ys, l_n);
    
    //Evaluation of spline
    for(int i = 0; i<N_eval+1; i++){
        double z = xs[0] + (xs[l_n-1]-xs[0])/N_eval*i;
        double spline_y = linterp(l_n, xs, ys, z);
        double integrated = linterp_integ(l_n, xs, ys, z);
        double g_interp = gsl_interp_eval(linear, xs, ys, z, NULL);
        double g_integ = gsl_interp_eval_integ(linear, xs, ys, xs[0], z, NULL);
        fprintf(loutput, "%g %g %g %g %g\n", z, spline_y, integrated, g_interp, g_integ);
    }
    gsl_interp_free(linear);
    fclose(linput);
    fclose(loutput);
    free(xs);
    free(ys);
    
    //B: QUADRATIC SPLINE
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
    
    //C: CUBIC SPLINE
    FILE *coutput = fopen("cubic.out", "w");
    //Define points from function
    int c_n = 20;
    xs = malloc(c_n*sizeof(double));
    ys = malloc(c_n*sizeof(double));
    
    fprintf(coutput, "#index 0 - points\n");
    for(int i = 0; i<c_n; i++){
        xs[i] = i;
        ys[i] = i + 10.*pow(sin(2.*i/M_PI), 5);
        fprintf(coutput, "%g %g\n", xs[i], ys[i]);
    }
    //Init cubic spline
    cspline *cs = cspline_init(c_n, xs, ys);

    fprintf(coutput, "\n\n#index 1 - interpolation\n");
    for(int i = 0; i<N_eval; i++){
        double z = xs[0] + (xs[c_n-1]-xs[0])/N_eval*i;
        double func = z + 10.*pow(sin(2.*z/M_PI), 5);
        fprintf(stderr, "test %d\n", i);
        double spline_y = cspline_eval(cs, z);
        
        fprintf(coutput, "%g %g %g\n", z, func, spline_y);
    }
    
    cspline_free(cs);
    fclose(coutput);
    free(xs);
    free(ys);
    return 0;
}