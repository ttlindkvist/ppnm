#include<stdio.h>
#include<math.h>
#include<limits.h>
#include<float.h>

int equals(double a, double b, double tau, double epsilon);

void name_digit(int i){
    switch(i){
        case 0:
            printf("zero\n");
            break;
        case 1:
            printf("one\n");
            break;
        case 2:
            printf("two\n");
            break;
        case 3:
            printf("three\n");
            break;
        case 4:
            printf("four\n");
            break;
        case 5:
            printf("five\n");
            break;
        case 6:
            printf("six\n");
            break;
        case 7:
            printf("seven\n");
            break;
        case 8:
            printf("eight\n");
            break;
        case 9:
            printf("nine\n");
            break;
        default:
            printf("not a digit\n");
            break;
    }
}

int main(){
    int i=1;
    while(i+1>i) {
        i++;
    }
    printf("My max int ( while ): %i\n",i);
    
    // for(i = 1; i<(i+1); i++){
    //     //Do nothing, the for loop does everything
    // }
    // printf("My max int ( for ): %i\n",i);

    // i = 0;
    // do{
    //     i++;
    // } while(i+1>i);
    // printf("My max int ( do while ): %i\n",i);
    printf("INT_MAX: %d\n", INT_MAX);

    i=0;
    while(i-1<i) {
        i--;
    }
    printf("My min int ( while ): %i\n",i);

    // for(i = 0; i>(i-1); i--){
    //     //Do nothing, the for loop does everything
    // }
    // printf("My min int ( for ): %i\n",i);

    // i = 0;
    // do{
    //     i--;
    // } while(i-1<i);
    // printf("My min int ( do while ): %i\n",i);
    printf("INT_MIN: %d\n", INT_MIN);


    //Machine epsilon
    double x=1;
    while(1+x!=1){x/=2;}
    x*=2;
    printf("Calculated double epsilon (while) %g\n", x);

    /*double y;
    for(y=1; 1+y!=1; y/=2){}
    y*=2;
    printf("Calculated double epsilon (for) %g\n", x);

    double z;
    for(y=1; 1+y!=1; y/=2){}
    y*=2;
    printf("Calculated double epsilon (for) %g\n", x);

    printf("DBL_EPSILON=%g\n", DBL_EPSILON);*/
    
    float fx=1;
    while(1+fx!=1){fx/=2;}
    fx*=2;
    printf("Calculated float epsilon (while) %g\n", fx);
    printf("FLT_EPSILON=%g\n", FLT_EPSILON);

    long double ldx = 1;
    while(1+ldx!=1){ldx/=2;}
    ldx*=2;
    printf("Calculated long double epsilon (while) %Lg\n", ldx);
    printf("LDBL_EPSILON=%Lg\n", LDBL_EPSILON);
    
    printf("Other methods (for and do while) commented out\n\n");


    //Harmonic series
    int max = INT_MAX;
    float sum_up_float = 1.0f;
    for(int i = 2; i<max; i++){
        sum_up_float += 1.0f/i;
    }
    printf("sum up float=%g\n", sum_up_float);

    float sum_down_float = 0.f;
    for(int i = max-1; i>0; i--){
        sum_down_float += 1.0f/i;
    }
    printf("sum down float=%g\n", sum_down_float);

    printf("Large discrepency since the float epsilon is rather large. \
    Summing up after a while, we are adding a small number to a large one, resulting in nothing done.\n\
    But summing down means the following addition becomes more on the order of the previous sum, making addition possible.\n\n");

    printf("This is the harmonic series, which does not converge.\n");
    double sum_up_double = 1.0;
    for(int i = 2; i<max; i++){
        sum_up_double += 1.0/i;
    }
    printf("sum up double=%g\n", sum_up_double);

    double sum_down_double = 0.;
    for(int i = max-1; i>0; i--){
        sum_down_double += 1.0/i;
    }
    printf("sum down double=%g\n", sum_down_double);


    double a = 100.3;
    double b = 102;

    printf("a=%g\nb=%g\nequal, tau=4 = %s\n", a, b, equals(a, b, 4, 0) == 1 ? "true" : "false");

    printf("name_digit(9)=");
    name_digit(9);

    return 0;
}