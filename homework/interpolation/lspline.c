#include<assert.h>
#include"spline.h"

double linterp(int n, double *x, double *y, double z){
    assert(n > 1 && z >= x[0] && z<=x[n-1]);
    int i =0, j=n-1; 
    while(j-i>1){ //Binary search
        int m = (i+j)/2;
        if(z>x[m]){
            i=m;
        } else {
            j=m;
        }
    }
    return y[i] + (y[i+1]-y[i])/(x[i+1]-x[i])*(z-x[i]);
}
double linterp_integ(int n, double *x, double *y, double z){
    double sum = 0;
    
    int i =0, j=n-1; 
    while(j-i>1){ //Binary search
        int m = (i+j)/2;
        if(z>x[m]){
            i=m;
        } else {
            j=m;
        }
    }
    //Geometric interpretation: Add up squares plus triangles under graphs
    for(int k = 0; k<i; k++){
        sum += (x[k+1]-x[k])*y[k];
        sum += 0.5 * (x[k+1]-x[k])*(y[k+1]-y[k]);
    }
    double spline_y = y[i] + (y[i+1]-y[i])/(x[i+1]-x[i])*(z-x[i]);
    sum += (z-x[i])*y[i];
    sum += 0.5 * (z-x[i])*(spline_y-y[i]);
    
    return sum;
}