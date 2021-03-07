#include<math.h>
#include<stdio.h>
#include<stdlib.h>

double ex(double x, int only_positive){
    if(x<0 && only_positive==1)     return 1/ex(-x, only_positive);
    if(fabs(x)>1./8)  return pow(ex(x/2, only_positive),2);
    return 1+x*(1+x/2*(1+x/3*(1+x/4*(1+x/5*(1+x/6*(1+x/7*(1+x/8*(1+x/9*(1+x/10)))))))));
}
double ordinary_sum_ex(double x){
    if(x<0) return 1/ordinary_sum_ex(-x);
    if(x>1./8)  return pow(ordinary_sum_ex(x/2), 2);
    double sum = 0;
    for(int i = 0; i<=10; i++){
        sum += pow(x, i)/tgamma(i+1);
    }
    // for(int i = 10; i>=0; i--){
    //     sum += pow(x, i)/tgamma(i+1);
    // }
    return sum;
}

int main(){
    FILE *output = fopen("out.dat", "w");
    
    int N = 10000;
    double xmin = -100, xmax = 100;
    for(int i = 0; i<N; i++){
        double x = xmin + (xmax-xmin)/N*i;
        double math_exp = exp(x);
        double our_exp = ex(x, 1);
        double our_exp_w_neg = ex(x, 0);
        double sum_exp = ordinary_sum_ex(x);
        fprintf(output, "%g\t%g\t%g\t%g\t%g\t%g\n", x, our_exp, math_exp, fabs(math_exp-our_exp)/math_exp, \
        fabs(our_exp-our_exp_w_neg)/our_exp, fabs(our_exp-sum_exp)/our_exp); 
    }
    fclose(output);
    return 0;
}