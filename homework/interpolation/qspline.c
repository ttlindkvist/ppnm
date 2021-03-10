#include"spline.h"
#include<stdlib.h>
#include<assert.h>

qspline *qspline_init(int n, double *x, double *y){
    qspline *s = (qspline*)malloc(sizeof(qspline));
    s->x = malloc(n*sizeof(double));
    s->y = malloc(n*sizeof(double));
    for(int i = 0; i<n; i++){
        s->x[i] = x[i];
        s->y[i] = y[i];
    }
    
    s->b = malloc((n-1)*sizeof(double));
    s->c = malloc((n-1)*sizeof(double));
    s->n = n;
    double p[n-1], h[n-1];
    for(int i = 0; i<n-1; i++){
        h[i] = x[i+1]-x[i];
        p[i] = (y[i+1]-y[i])/h[i];
    }
    //Compute c coefficients recusively, first forward with guess c0=0
    s->c[0] = 0;
    for(int i = 0; i<n-2; i++){
        s->c[i+1] = (p[i+1]-p[i]-s->c[i]*h[i])/h[i+1];
    }
    //Now backwards with half of c_n-1 from before 
    s->c[n-2] /= 2.;
    for(int i = n-3; i>=0; i--){
        s->c[i] = (p[i+1]-p[i]-s->c[i+1]*h[i+1])/h[i];
    }
    
    //Calculate b coefficients
    for(int i = 0; i<n-1; i++){
        s->b[i] = p[i]-s->c[i]*h[i];
    }
    return s;
}
double qspline_eval(qspline *s, double z){
    assert(z>= s->x[0] && z<=s->x[s->n-1]);
    //Binary search
    int i =0, j=s->n-1; 
    while(j-i>1){ //Binary search
        int m = (i+j)/2;
        if(z>s->x[m]){
            i=m;
        } else {
            j=m;
        }
    }
    double h = z-s->x[i];
    return s->y[i] + h*(s->b[i]+h*s->c[i]);
}
double qspline_deriv(qspline *s, double z){
    assert(z>= s->x[0] && z<=s->x[s->n-1]);
    //Binary search
    int i =0, j=s->n-1; 
    while(j-i>1){ //Binary search
        int m = (i+j)/2;
        if(z>s->x[m]){
            i=m;
        } else {
            j=m;
        }
    }
    double h = z-s->x[i];
    return s->b[i]+2*h*s->c[i];
}
double qspline_integ(qspline *s, double z){
    assert(z>= s->x[0] && z<=s->x[s->n-1]);
    //Binary search
    int i =0, j=s->n-1; 
    while(j-i>1){ //Binary search
        int m = (i+j)/2;
        if(z>s->x[m]){
            i=m;
        } else {
            j=m;
        }
    }
    double sum = 0;
    for(int k = 0; k<i; k++){
        double delta_x = s->x[k+1]-s->x[k];
        sum += s->y[k]*delta_x;
        sum += 0.5 * s->b[k]*delta_x*delta_x;
        sum += 1./3 * s->c[k]*delta_x*delta_x*delta_x;
    }
    
    double h = z-s->x[i];
    sum += h*s->y[i];
    sum += 0.5 * s->b[i]*h*h;
    sum += 1./3 * s->c[i]*h*h*h;
    return sum;
}
void qspline_free(qspline *s){
    free(s->b);
    free(s->c);
    free(s);    
}