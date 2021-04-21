#include<iostream>
#include<vector>
#include<array>
#include<cmath>
#include"montecarlo.hpp"

double f_inv_cosines(std::vector<double> &x){
    return 1./(1.-cos(x[0])*cos(x[1])*cos(x[2]));
}

void test_function(double f(std::vector<double> &), std::vector<double> &a, std::vector<double> &b, int N, double exact){
    
    std::pair<double, double> I1 = plainmc(f, a, b, N);
    std::cout << I1.first << std::endl;
    std::cout << "Should be: " << exact << std::endl;
}

int main(){
    std::vector<double> a{0, 0, 0};
    std::vector<double> b{M_PI, M_PI, M_PI};
    
    test_function(f_inv_cosines, a, b, 10000000, 1.3932039296856*M_PI*M_PI*M_PI);

    
    return 0;
}