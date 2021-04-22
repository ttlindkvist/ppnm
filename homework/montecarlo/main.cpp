#include<iostream>
#include<iomanip>
#include<vector>
#include<array>
#include<cmath>
#include"montecarlo.hpp"

double f_inv_cosines(double *x){
    return 1./(1.-cos(x[0])*cos(x[1])*cos(x[2]));
}
double vol_sphere(double *x){
    double r2 = x[0]*x[0]+x[1]*x[1]+x[2]*x[2];
    return r2 <= 1 ? 1 : 0; 
}
double area_circle(double *x){
    return x[0]*x[0]+x[1]*x[1] <= 1 ? 1 : 0;
}

void test_function_plain(double f(double *x), std::vector<double> &a, std::vector<double> &b, int N, double exact){
    std::pair<double, double> I1 = plainmc(f, a, b, N);
    std::cout << "Calculated to:   " << I1.first  << std::endl;
    std::cout << "Exact result:    " << exact     << std::endl;
    std::cout << "Estimated error: " << I1.second << std::endl;
    std::cout << "Actual error:    " << fabs(I1.first-exact) << std::endl << std::endl;
}
void test_function_quasi(double f(double *x), std::vector<double> &a, std::vector<double> &b, int N, double exact){
    std::pair<double, double> I1 = quasi_mc(f, a, b, N);
    std::cout << "Calculated to:   " << I1.first  << std::endl;
    std::cout << "Exact result:    " << exact     << std::endl;
    std::cout << "Estimated error: " << I1.second << std::endl;
    std::cout << "Actual error:    " << fabs(I1.first-exact) << std::endl << std::endl;
}
void partA(int N){
    std::vector<double> a{-1, -1};
    std::vector<double> b{1, 1};
    std::cout << "Calculating area of circle\nN = " << N << std::endl;
    test_function_plain(area_circle, a, b, N, M_PI);

    a.push_back(-1);
    b.push_back(1);

    std::cout << "Calculating volume of sphere with radius 1\nN = " << N << std::endl;
    test_function_plain(vol_sphere, a, b, N, 4.*M_PI/3.);

    a = std::vector<double>{0, 0, 0};
    b = std::vector<double>{M_PI, M_PI, M_PI};

    std::cout << "Calculating int (1-cos(x)cos(y)cos(z))^-1 with x,y,z=(0, pi)\nN = " << N << std::endl;
    test_function_plain(f_inv_cosines, a, b, N, 1.3932039296856*M_PI*M_PI*M_PI);

}
void partB(int N){
    std::vector<double> a{-1, -1};
    std::vector<double> b{1, 1};
    
    std::cout << "Calculating area of circle\nN = " << N << std::endl;
    test_function_quasi(area_circle, a, b, N, M_PI);

    a.push_back(-1);
    b.push_back(1);

    std::cout << "Calculating volume of sphere with radius 1\nN = " << N << std::endl;
    test_function_quasi(vol_sphere, a, b, N, 4.*M_PI/3.);

    a = std::vector<double>{0, 0, 0};
    b = std::vector<double>{M_PI, M_PI, M_PI};

    std::cout << "Calculating int (1-cos(x)cos(y)cos(z))^-1 with x,y,z=(0, pi)\nN = " << N << std::endl;
    test_function_quasi(f_inv_cosines, a, b, N, 1.3932039296856*M_PI*M_PI*M_PI);
}
void compareError(int Na, int Nb, int step, FILE *out){

    std::vector<double> a{-1, -1};
    std::vector<double> b{1, 1};
    for(int i = Na; i<Nb; i+=step){
        std::pair<double,double> plain = plainmc(area_circle, a, b, i);
        std::pair<double,double> quasi = quasi_mc(area_circle, a, b, i);
        
        fprintf(out, "%d %g %g %g %g\n", i, fabs(plain.first-M_PI), fabs(quasi.first-M_PI), plain.second, quasi.second);
    }
}
int main(){
    std::cout << std::fixed;
    std::cout << std::setprecision(15);
    srand(42);
    
    std::cout << "Monte Carlo integration\n\n----------- Part A - plain monte carlo integration\n\n";
    partA(1e6);
    
    
    std::cout << "\n----------- Part B ---------\nSame functions but now this low-discrepancy (Halton) sequences\n\n";
    partB(1e6);
   
    std::cout << "Error comparison on volume of sphere\n";
    FILE *output = fopen("error.out", "w");
    compareError(100, 2e5, 1e3, output);

    fclose(output);
    return 0;
}