#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<assert.h>
#include<omp.h>
#include<gsl/gsl_qrng.h>

typedef struct {int N; unsigned int seed; int N_in;} monte_carlo_data;

void *monte_carlo_pi(void * a){
    monte_carlo_data *args = (monte_carlo_data*)a;
    int N = args->N;
    unsigned int seed = args->seed;

    int N_in = 0;
    for(int i = 0; i<N; i++){
        double x = (double)rand_r(&seed)/RAND_MAX;
        double y = (double)rand_r(&seed)/RAND_MAX;
        if(x*x+y*y < 1.0){
            N_in++;
        }
    }
    args->N_in = N_in;
    return NULL;
}
void *monte_carlo_pi_gsl(void * a){
    monte_carlo_data *args = (monte_carlo_data*)a;
    int N = args->N;
    int N_in = 0;

    gsl_qrng * q = gsl_qrng_alloc (gsl_qrng_sobol, 2);

    for(int i = 0; i<N; i++){
        double v[2];
        gsl_qrng_get(q, v);
        if(v[0]*v[0]+v[1]*v[1] < 1.0){
            N_in++;
        }
    }
    args->N_in = N_in;
    return NULL;
}

int main(int argc, char **argv){
    assert(argc == 2);
    int nthreads = atoi(argv[1]);
    int reps = 2000;
    unsigned int curr_time = time(NULL);
    double N_in_values[reps];
    double N_in_values_gsl[reps];

    FILE *output = fopen("convergence.dat", "w");
    #pragma omp parallel for
    for(int t = 0; t<nthreads; t++){
        for(int i = 1; i<reps/nthreads; i++){
            int N = pow(10, 6.0/reps*(i*nthreads+t));
            monte_carlo_data dat = {.N = N, .seed = curr_time+t+i};
            monte_carlo_data dat_gsl = {.N = N, .seed = curr_time+t+i};
            monte_carlo_pi((void*)&dat);
            monte_carlo_pi_gsl((void*)&dat_gsl);
            N_in_values[i*nthreads + t] = dat.N_in;
            N_in_values_gsl[i*nthreads + t] = dat_gsl.N_in;
        }
    }
    for(int t = 0; t<nthreads; t++){
        for(int i = 1; i<=reps/nthreads; i++){
            int N = pow(10, 6.0/reps*(i*nthreads+t));
            fprintf(output, "%g\t%g\t%g\t%g\n", (double)N, (4.0*N_in_values[i*nthreads+t]/N-M_PI)/M_PI,(4.0*N_in_values_gsl[i*nthreads+t]/N-M_PI)/M_PI, 1.0/sqrt(N));
        }
    }
    
    fclose(output);
    return 0;
}