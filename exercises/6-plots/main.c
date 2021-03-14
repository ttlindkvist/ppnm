#include<stdio.h>
#include<math.h>
#include<complex.h>
#include<gsl/gsl_sf_erf.h>
#include<gsl/gsl_sf_gamma.h>

double Erf(double x){
    /// single precision error function (Abramowitz and Stegun, from Wikipedia)
    if(x<0) return -Erf(-x);
    double a[]={0.254829592,-0.284496736,1.421413741,-1.453152027,1.061405429};
    double t=1/(1+0.3275911*x);
    double sum=t*(a[0]+t*(a[1]+t*(a[2]+t*(a[3]+t*a[4]))));/* the right thing */
    return 1-sum*exp(-x*x);
}
double Gamma(double x){
    ///single precision gamma function (Gergo Nemes, from Wikipedia)
    if(x<0)return M_PI/sin(M_PI*x)/Gamma(1-x);
    if(x<9)return Gamma(x+1)/x;
    double lnGamma=x*log(x+1/(12*x-1/x/10))-x+log(2*M_PI/x)/2;
    return exp(lnGamma);
}
double LGamma(double x){
    ///single precision gamma function (Gergo Nemes, from Wikipedia);
    if(x<9)return LGamma(x+1)-log(x);
    double lnGamma=x*log(x+1/(12*x-1/x/10))-x+log(2*M_PI/x)/2;
    return lnGamma;
}

complex CGamma(complex z){
    if(creal(z)<0)return M_PI/csin(M_PI*z)/CGamma(1-z);
    if(creal(z)<9)return CGamma(z+1)/z;
    double lnGamma=z*clog(z+1/(12*z-1/z/10))-z+clog(2*M_PI/z)/2;
    return cexp(lnGamma);
}

int main(){
    FILE *erf_out = fopen("erf.dat", "w");
    FILE *gamma_out = fopen("gamma.dat", "w");
    FILE *lgamma_out = fopen("lgamma.dat", "w");
    FILE *cgamma_out = fopen("cgamma.dat", "w"); 

    double xmin=-2, xmax=2;
    for(double x=xmin; x<=xmax; x+=1.0/16){
        fprintf(erf_out, "%10g %10g %10g %10g\n", x, erf(x), gsl_sf_erf(x), Erf(x));
    }

    int N = 10000;
    double eps = 1.0/(N*10);
    for(int i=0; i<N; i++){
        double x = -5.0 + 10.0*i/N + eps;
        fprintf(gamma_out, "%10g %10g %10g %10g\n", x, tgamma(x), gsl_sf_gamma(x), Gamma(x));
    }

    for(int i=0; i<N; i++){
        double x = 10.0*i/N + eps;
        fprintf(lgamma_out, "%10g %10g %10g %10g\n", x, lgamma(x), gsl_sf_lngamma(x), LGamma(x));
    }

    N = 200;
    for(int i=0; i<N; i++){
        fprintf(cgamma_out, "\n");
        for(int j=0; j<N; j++){
            double x = -4.0 + 8.0*i/N + eps;
            double y = -4.0 + 8.0*j/N + eps;
            complex g = CGamma(x+I*y);
            fprintf(cgamma_out, "%10g %10g %10g\n", x, y, cabs(g));
        }
        
    }

    fclose(erf_out);
    fclose(gamma_out);
    fclose(lgamma_out);
    fclose(cgamma_out);
    return 0;
}