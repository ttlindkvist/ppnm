#include<assert.h>
#include<stdlib.h>
#include"spline.h"

cspline *cspline_init(int n, double *x, double *y){
    //Allocate memory
    cspline *s = (cspline*)malloc(sizeof(cspline));
    s->n = n;
    s->x = malloc(n*sizeof(double));
    s->y = malloc(n*sizeof(double));
    for(int i = 0; i<n; i++){
        s->x[i] = x[i];
        s->y[i] = y[i];
    }
    s->b = malloc(n*sizeof(double));
    s->c = malloc((n-1)*sizeof(double));
    s->d = malloc((n-1)*sizeof(double));
    
    //Define usefull values h=delta x and p = linear slope coefficient "lineær hældningskoefficient" 
    double h[n-1], p[n-1];
    for(int i = 0; i<n-1; i++){
        h[i] = x[i+1]-x[i];
        assert(h[i] > 0);
        p[i] = (y[i+1]-y[i])/h[i];
    }
    //Building the tri-diagonal system
    double D[n], Q[n-1], B[n];
    D[0] = 2; D[n-1] = 2; Q[0] = 1; B[0] = 3*p[0]; B[n-1] = 3*p[n-2];
    for(int i = 0; i<n-2; i++){
        double frac_h = h[i]/h[i+1];
        D[i+1] = 2*frac_h + 2;
        Q[i+1] = frac_h;
        B[i+1] = 3*(p[i]+p[i+1]*frac_h);
    }
    //Gauss-elimination to upper-triangular matrix
    for(int i = 1; i<n; i++){
        D[i] -= Q[i-1]/D[i-1];
        B[i] -= B[i-1]/D[i-1];
    }
    //Back-substitution
    s->b[n-1]=B[n-1]/D[n-1];
    for(int i = n-2; i>=0; i--){
        s->b[i] = (B[i]-Q[i]*s->b[i+1])/D[i];
    }
    //Calculate c and d coefficients
    for(int i = 0; i<n-1; i++){
        s->c[i] = (-2*s->b[i] - s->b[i+1]+3*p[i])/h[i];
        s->d[i] = (s->b[i] + s->b[i+1] - 2*p[i])/(h[i]*h[i]);
    }
    
    return s;
}
void cspline_free(cspline *s){
    free(s->b);
    free(s->c);
    free(s->d);
    free(s->x);
    free(s->y);
    free(s);
}
double cspline_eval(cspline *s, double z){
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
    return s->y[i] + s->b[i]*h + s->c[i]*h*h + s->d[i]*h*h*h;
}
double cspline_deriv(cspline *s, double z){
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
    return s->b[i] + 2*s->c[i]*h + 3*s->d[i]*h*h;
}
double cspline_integ(cspline *s, double z){
    return 0;
}