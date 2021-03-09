#include<stdio.h>
#include<math.h>
#include<assert.h>
#include<stdlib.h>

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
    //Add up squares plus triangles under graphs
    for(int k = 0; k<i; k++){
        sum += (x[k+1]-x[k])*y[k];
        sum += 0.5 * (x[k+1]-x[k])*(y[k+1]-y[k]);
    }
    double spline_y = y[i] + (y[i+1]-y[i])/(x[i+1]-x[i])*(z-x[i]);
    sum += (z-x[i])*y[i];
    sum += 0.5 * (z-x[i])*(spline_y-y[i]);
    return sum;
}

int main(){
    FILE *output = fopen("linear.out", "w");
    FILE *input = fopen("points.in", "r");
    
    int npoints;
    int ret = fscanf(input, "%d", &npoints);
    assert(ret>0);
    assert(npoints > 1);
    
    double *xs = calloc(npoints, sizeof(double));
    double *ys = calloc(npoints, sizeof(double));    
    double x, y;
    int i = 0;
    while(fscanf(input, "%lg %lg", &x, &y) > 0){
        assert(i<npoints);
        fprintf(stderr, "%g %g\n", x, y);
        xs[i] = x;
        ys[i] = y;
        i++;
    }
    
    int N_eval = 200;
    for(int i = 0; i<N_eval+1; i++){
        double z = xs[0] + (xs[npoints-1]-xs[0])/N_eval*i;
        double liny = linterp(npoints, xs, ys, z);
        double integrated = linterp_integ(npoints, xs, ys, z);
        fprintf(output, "%g %g %g\n", z, liny, integrated);
    }
    fclose(output);
    fclose(input);
    free(xs);
    free(ys);
    return 0;
}