#include<math.h>
#include<stdio.h>
#include<stdlib.h>

double ex(double x, int only_positive){
    if(x<0 && only_positive==1)     return 1/ex(-x, only_positive);
    if(fabs(x)>1./8)  return pow(ex(x/2, only_positive),2);
    return 1+x*(1+x/2*(1+x/3*(1+x/4*(1+x/5*(1+x/6*(1+x/7*(1+x/8*(1+x/9*(1+x/10)))))))));
}
double ordinary_sum_ex(double x, int forwards){
    if(x<0) return 1/ordinary_sum_ex(-x, forwards);
    if(x>1./8)  return pow(ordinary_sum_ex(x/2, forwards), 2);
    double sum = 0;
    if(forwards == 1){
        for(int i = 0; i<=10; i++){
            sum += pow(x, i)/tgamma(i+1);
        }
    } else {
        for(int i = 10; i>=0; i--){
            sum += pow(x, i)/tgamma(i+1);
        }
    }
    return sum;
}

int main(){
    FILE *output = fopen("out.dat", "w");
    FILE *output_large_N = fopen("largeN.dat", "w");
    int N = 200;
    double xmin = -100, xmax = 100;
    for(int i = 0; i<N; i++){
        double x = xmin + (xmax-xmin)/N*i;
        double math_exp = exp(x);
        double our_exp = ex(x, 1);
        double our_exp_w_neg = ex(x, 0);
        double sum_exp_backwards = ordinary_sum_ex(x, 0);
        double sum_exp_forwards = ordinary_sum_ex(x, 1);
        fprintf(output, "%g\t%g\t%g\t%g\t%g\t%g\t%g\n", x, our_exp, math_exp, fabs(math_exp-our_exp)/math_exp, \
        fabs(math_exp-our_exp_w_neg)/math_exp, fabs(math_exp-sum_exp_backwards)/math_exp, fabs(math_exp-sum_exp_forwards)/math_exp);
    }
    int largeN = 10000;
    for(int i = 0; i<largeN; i++){
        double x = xmin + (xmax-xmin)/largeN*i;
        double math_exp = exp(x);
        double our_exp = ex(x, 1);
        fprintf(output_large_N, "%g\t%g\n", x, fabs(math_exp-our_exp)/math_exp);
    }
    
    fclose(output);
    fclose(output_large_N);
    return 0;
}