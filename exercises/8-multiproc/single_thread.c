#include<pthread.h>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<assert.h>

typedef struct {int N; unsigned int *seed; int N_in;} monte_carlo_data;

void *monte_carlo_pi(void * a){
    monte_carlo_data *args = (monte_carlo_data*)a;
    int N = args->N;
    unsigned int *seed = args->seed;

    int N_in = 0;
    for(int i = 0; i<N; i++){
        double x = 1.0*rand_r(seed)/RAND_MAX;
        double y = 1.0*rand_r(seed)/RAND_MAX;
        if(x*x+y*y < 1.0){
            N_in++;
        }
    }
    args->N_in = N_in;
    return NULL;
}


int main(int argc, char **argv){
    assert(argc>=2);
    int N = atoi(argv[1]); unsigned int seed = 1;  
    monte_carlo_data dat = {.N = N, .seed = &seed};
    monte_carlo_pi((void*)&dat);
    double final_pi = 4.0*dat.N_in/(N);
    printf("\nSingle thread with N=%g\nMC value of pi = %.15g\n", (double)N, final_pi);
    
    return 0;
}