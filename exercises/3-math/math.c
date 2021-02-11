#include<stdio.h>
#include<math.h>
#include<complex.h>

int main(){
	printf("gamma(5)=%g\n", tgamma(5));
	printf("J_1(0.5)=%.10g\n", j1(0.5));
	complex z = csqrt(-2);
	printf("sqrt(-2)=%g + i * %g\n", creal(z), cimag(z));

	z = cpow(M_E, I*M_PI);
	printf("e^(i*pi) = %g + i * %g\n", creal(z), cimag(z));

	z = cpow(M_E, I);
	printf("e^i = %g + i * %g\n", creal(z), cimag(z));

	z = cpow(I, M_E);
	printf("i^e = %g + i * %g\n", creal(z), cimag(z));

	z = cpow(I, I);
	printf("i^i = %g + i * %g\n", creal(z), cimag(z));

	float x_float = 1.f/9;
	double x_double = 1./9;
	long double x_long_double = 1.L/9;

	printf("Checking of significant digits of float, double and long double\n");
	printf("%.25g\n%.25lg\n%.25Lg\n", x_float, x_double, x_long_double);
return 0;
}
