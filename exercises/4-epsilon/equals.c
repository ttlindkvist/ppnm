#include <math.h>

//Returns 1 if a and b are equal with absolure precision 'tau'
// or with relative precision 'epsilon'
int equals(double a, double b, double tau, double epsilon){
    if(fabs(a-b)<tau){
        return 1;
    }
    if(fabs(a-b)/(fabs(a)+fabs(b)) < epsilon/2){
        return 1;
    }
    return 0;
}