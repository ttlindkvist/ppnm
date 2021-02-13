#include "komplex.h"
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<complex.h>
#define TAU 1e-6
#define EPS 1e-6

#define KOMPLEX(z) komplex_new(creal(z), cimag(z))

int main(){
	srand(time(NULL));
	komplex a = {rand() % 100 - 50, rand() % 100 - 50};
	komplex b = {rand() % 100 - 50, rand() % 100 - 50};
	complex A = a.re + a.im * I;
	complex B = b.re + b.im * I;

	komplex_print("a =", a);
	komplex_print("b =", b);
	komplex_print("a+b =", komplex_add(a, b));
	komplex_print("a*b =", komplex_mul(a, b));
	komplex_print("a/b =", komplex_div(a, b));

	komplex_print("cos(a)=", komplex_cos(a));
	if(komplex_equal(KOMPLEX(ccos(A)), komplex_cos(a),TAU,EPS)){
		printf("komplex cos passed \n");
	} else {
		printf("komplex cos failed \n");
	}

	komplex_print("sin(b)=", komplex_sin(b));
	if(komplex_equal(KOMPLEX(csin(B)), komplex_sin(b),TAU,EPS)){
		printf("komplex sin passed \n");
	} else {
		printf("komplex sin failed \n");
	}

	komplex_print("sqrt(a)=", komplex_sqrt(a));
	if(komplex_equal(KOMPLEX(csqrt(A)), komplex_sqrt(a),TAU,EPS)){
		printf("komplex sqrt passed \n");
	} else {
		printf("komplex sqrt failed \n");
	}
	return 0;
}
