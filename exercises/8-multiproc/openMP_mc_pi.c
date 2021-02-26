#include<stdio.h>
#include<math.h>
#include<omp.h>
#include<stdlib.h>
#include<time.h>
#include<assert.h>

typedef struct {int N; unsigned int seed; int N_in;} monte_carlo_data;

void *monte_carlo_pi(void * a){
    monte_carlo_data *args = (monte_carlo_data*)a;
    int N = args->N;
    unsigned int seed = args->seed;

    int N_in = 0;
    for(int i = 0; i<N; i++){
        double x = 1.0*rand_r(&seed)/RAND_MAX;
        double y = 1.0*rand_r(&seed)/RAND_MAX;
        if(x*x+y*y < 1.0){
            N_in++;
        }
    }
    args->N_in = N_in;
    return NULL;
}

int main(int argc, char **argv){
    assert(argc > 2);
    int nthreads = atoi(argv[1]);
    int N = atoi(argv[2]);
    int N_in[nthreads];
    int t = time(NULL);
    unsigned int seeds[nthreads];

    printf("\nopenmp with %d threads and total N = %g\n", nthreads, (double)N*nthreads);


    #pragma omp parallel for
    for(int i = 0; i<nthreads; i++){
        seeds[i] = t+2*i;
        monte_carlo_data dat = {.N = N, .seed = seeds[i]};
        monte_carlo_pi((void*)&dat);
        N_in[i] = dat.N_in;
    }

    int total_N_in = 0;
    for(int i = 0; i<nthreads; i++){
        total_N_in += N_in[i];
    }
    double final_pi = 4.0*total_N_in/(N*nthreads);
    printf("MC value of pi = %.15g\n\n", final_pi);
    return 0;
}