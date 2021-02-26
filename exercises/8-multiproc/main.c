#include<pthread.h>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

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
    int nthreads = 4;
    int N = 1e7;

    int t = time(NULL);
    if(argc > 1){
        nthreads = atoi(argv[1]);
        N = atoi(argv[2]);
    }
    pthread_t threads[nthreads];
    pthread_attr_t *attributes = NULL;
    int flags[nthreads];
    unsigned int seeds[nthreads];
    monte_carlo_data dat[nthreads];

    for(int i = 0; i<nthreads; i++){
        seeds[i] = t+i;
        dat[i].N = N;
        dat[i].seed = &seeds[i];
        flags[i] = pthread_create(&threads[i], attributes, monte_carlo_pi, (void*)&dat[i]);
    }
    for(int i = 0; i<nthreads; i++){
        flags[i] = pthread_join(threads[i], NULL);
        printf("Thread nr %d returns %g\n", i, 4.0*dat[i].N_in/dat[i].N);
    }
    
    return 0;
}