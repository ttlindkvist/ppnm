#include<stdio.h>
#include<math.h>
#include<omp.h>
#include<stdlib.h>
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
    assert(argc > 1);
    int nthreads = 4;
    int N = atoi(argv[1]);
    int total_N_in = 0;
    int t = time(NULL);
    unsigned int seeds[nthreads];

    printf("\nopenmp with %d threads and total N = %g\n", nthreads, (double)N);


    #pragma omp parallel sections
    {
        #pragma omp section
        {
            seeds[0] = t;
            monte_carlo_data dat = {.N = N/4, .seed = &seeds[0]};
            monte_carlo_pi((void*)&dat);
            total_N_in += dat.N_in;
        }
        #pragma omp section
        {
            seeds[1] = t+1;
            monte_carlo_data dat = {.N = N/4, .seed = &seeds[1]};
            monte_carlo_pi((void*)&dat);
            total_N_in += dat.N_in;

        }
        #pragma omp section
        {
            seeds[2] = t+2;
            monte_carlo_data dat = {.N = N/4, .seed = &seeds[2]};
            monte_carlo_pi((void*)&dat);
            total_N_in += dat.N_in;

        }
        #pragma omp section
        {
            seeds[3] = t+3;
            monte_carlo_data dat = {.N = N/4, .seed = &seeds[3]};
            monte_carlo_pi((void*)&dat);
            total_N_in += dat.N_in;

        }
    }
    double final_pi = 4.0*total_N_in/N;
    printf("Avg value of pi = %.15g\n\n", final_pi);
    return 0;
}