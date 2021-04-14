#include "integration.hpp"
#include<iostream>
#include<cmath>
#include<functional>

void test_function(std::function<std::pair<double,int>(double)> f){
    std::pair<double,int> f_call = f(10); 
    f_call = f(20);
    std::cout << f_call.first << " " << f_call.second << std::endl;
}

int main(){

    double eps = 1e-4, abs = 1e-4;

    printf("\n------ PART A ------ \nAll integrals calculated with recursive adaptive integrator\n");
    printf("Absolute precision %g - Relative precision %g\n", abs, eps);

    int fcalls = 0;
    double I = adapt_quad24([](double x){return sqrt(x);}, 0, 1, abs, eps);
    printf("int sqrt(x) from 0 to 1\nCalculated to \t%.16f\nShould be \t%.16f\nWith error\t%.16f\nwith %d function evaluations\n\n", \
            I, 2./3, fabs(2./3-I), fcalls);

    return 0;
}