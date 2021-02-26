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
    if(argc > 2){
        nthreads = atoi(argv[1]);
        N = atoi(argv[2]);
    }
    printf("\npthreads with %d threads and total N = %g\n", nthreads, (double)N*nthreads);
    pthread_t threads[nthreads-1];
    
    unsigned int seeds[nthreads];
    monte_carlo_data dat[nthreads];

    for(int i = 0; i<nthreads-1; i++){
        seeds[i] = t+i;
        dat[i].N = N;
        dat[i].seed = &seeds[i];
        pthread_create(&threads[i], NULL, monte_carlo_pi, (void*)&dat[i]);
    }

    //Utilize main thread
    seeds[nthreads-1] = t+nthreads-1;
    dat[nthreads-1].N = N;
    dat[nthreads-1].seed = &seeds[nthreads-1];
    monte_carlo_pi((void*)&dat[nthreads-1]);
    
    double total_N_in = dat[nthreads-1].N_in;
    for(int i = 0; i<nthreads-1; i++){
        pthread_join(threads[i], NULL);
        total_N_in += dat[i].N_in;
    }

    double final_pi = 4.0*total_N_in/(N*nthreads);
    printf("Avg value of pi = %.15g\n\n", final_pi);
    
    return 0;
}