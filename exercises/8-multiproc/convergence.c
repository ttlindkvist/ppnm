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
        double x = (double)rand_r(seed)/RAND_MAX;
        double y = (double)rand_r(seed)/RAND_MAX;
        if(x*x+y*y < 1.0){
            N_in++;
        }
    }
    args->N_in = N_in;
    return NULL;
}


int main(int argc, char **argv){
    int reps = 1000;
    unsigned int seed = time(NULL);

    FILE *output = fopen("convergence.dat", "w");

    for(int i = 1; i<=reps; i++){
        int N = pow(10, 6.0/reps*i);
        monte_carlo_data dat = {.N = N, .seed = &seed};
        monte_carlo_pi((void*)&dat);
        double final_pi = 4.0*dat.N_in/N;
        fprintf(output, "%g\t%g\t%g\t%g\n", (double)N, final_pi, (final_pi-M_PI)/M_PI, 1.0/sqrt(N));
    }
    fclose(output);
    return 0;
}