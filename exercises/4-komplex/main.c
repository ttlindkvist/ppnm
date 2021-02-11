#include "komplex.h"
#include<stdio.h>

int main(){
	komplex a = {1, 2};
	komplex b = {5, 2};
	komplex_print("a =", a);
	komplex_print("b =", b);
	komplex_print("a+b =", komplex_add(a, b));
	komplex_print("a*b =", komplex_mul(a, b));
	komplex_print("a/b =", komplex_div(a, b));
	return 0;
}
