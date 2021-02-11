#include "komplex.h"
#include <stdio.h>

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

// int komplex_equal(komplex a, komplex b, double acc, double eps){
//     komplex diff = komplex_sub(a, b);

//     if(diff.re < acc && diff.im < acc && ){
//     }
// }