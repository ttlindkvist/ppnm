#include "ann.hpp"
#include<cstdio>
#include<cmath>
#include<vector>

double gauss_wavelet(double x){
    return x*exp(-x*x);
}

int main(){
    std::vector<double> xs;
    std::vector<double> ys;
    double a = -1;
    double b = 1;
    double points = 20;
    for(int i = 0; i<points; i++){
        double x = a + (b-a)*i/points;
        double y = cos(5*x-1)*exp(-x*x);
        xs.push_back(x);
        ys.push_back(y);
    }
    
    gsl_vector_view x_view = gsl_vector_view_array(xs.data(), xs.size());
    gsl_vector_view y_view = gsl_vector_view_array(ys.data(), ys.size());
    
    ANN ann(4, gauss_wavelet);
    ann.train(&x_view.vector, &y_view.vector);
    
    
    FILE *out = fopen("interp.out", "w");
    int N = 100;
    fprintf(out, "#index 0 -- points\n");
    for(int i = 0; i<xs.size(); i++){
        fprintf(out, "%g %g\n", xs[i], ys[i]);
    }
    
    fprintf(out, "\n\n#index 1 -- interpolated\n");
    for(int i = 0; i<N; i++){
        double x = xs[0] + (xs[xs.size()-1]-xs[0])*1.0*i/N;
        fprintf(out, "%g %g\n", x, ann.response(x));
    }
    fclose(out);
    return 0;
}