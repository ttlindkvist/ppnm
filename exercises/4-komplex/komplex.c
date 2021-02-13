#include "komplex.h"
#include <stdio.h>
#include <math.h>

void komplex_print(const char *s, komplex a){
    printf("%s (%g, %g)\n", s, a.re, a.im);
}

void komplex_set(komplex *z, double x, double y){
    z->re = x;
    z->im = y;
}

komplex komplex_new(double x, double y){
    komplex z = {x, y};
    return z;
}

komplex komplex_add(komplex a, komplex b){
    komplex z = {a.re + b.re, a.im + b.im};
    return z;
}

komplex komplex_sub(komplex a, komplex b){
    komplex z = {a.re - b.re, a.im - b.im};
    return z;
}

int komplex_equal(komplex a, komplex b, double acc, double eps){
    komplex diff = komplex_sub(a, b);

    if(fabs(diff.re) < acc && fabs(diff.im) < acc && (fabs(diff.re)/(fabs(a.re) + fabs(b.re)) < eps/2 && \
        fabs(diff.im)/(fabs(a.im) + fabs(b.im)) < eps/2)){
        return 1;  
    }
    return 0;
}

komplex komplex_mul(komplex a, komplex b){
    komplex z = {a.re*b.re - a.im*b.im, a.re*b.im + a.im*b.re};
    return z; 
}

komplex komplex_div(komplex a, komplex b){
    komplex z = komplex_mul(a, komplex_conj(b));
    double abs2 = b.re*b.re + b.im*b.im;
    z.re /= abs2;
    z.im /= abs2;
    return z;
}

komplex komplex_conj(komplex a){
    komplex z = {a.re, -a.im};
    return z;
}

double komplex_abs(komplex z){
    return sqrt(komplex_mul(z, komplex_conj(z)).re);
}

komplex komplex_exp(komplex z){
    //exp(a+ib) = exp(a)*exp(ib) = exp(a)*(cos(b)+isin(b))
    double real_exp = exp(z.re);
    komplex r = {real_exp*cos(z.im), real_exp*sin(z.im)};
    return r;
}

komplex komplex_sin(komplex z){
    komplex im = {0, 1};
    komplex neg_im = {0, -1};
    komplex half_im = {0, -1.0/2.0};
    return komplex_mul(half_im, komplex_sub(komplex_exp(komplex_mul(im, z)), komplex_exp(komplex_mul(neg_im, z))));
}
komplex komplex_cos(komplex z){
    komplex im = {0, 1};
    komplex neg_im = {0, -1};
    komplex half = {1.0/2.0, 0};
    return komplex_mul(half, komplex_add(komplex_exp(komplex_mul(im, z)), komplex_exp(komplex_mul(neg_im, z))));
}

komplex komplex_sqrt(komplex z){
    double x = z.re;
    double y = z.im;
    int im_sign = y > 0 ? 1 : -1;
    komplex r = {sqrt((sqrt(x*x+y*y)+x)/2), im_sign*sqrt((sqrt(x*x+y*y)-x)/2)};
    return r;
}