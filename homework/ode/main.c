#include<stdio.h>
#include<assert.h>
#include<gsl/gsl_blas.h>
#include"ode.h"
#include"matrix.h"

void func_shm(double t, gsl_vector *y, gsl_vector *dydt){
    // assert(y->size == dydt->size && y->size == 2);
    //u''=-u
    //y0 = u
    //y1 = u'
    //
    //y0' = y1
    //y1' = -y0
    
    gsl_vector_set(dydt, 0, gsl_vector_get(y, 1));
    gsl_vector_set(dydt, 1, -gsl_vector_get(y, 0));
}
double Tc = 5.;    /* Time between contacts (days)*/
double Tr = 14.; /* Recovery time */
void func_sir(double t, gsl_vector *y, gsl_vector *dydt){
    // assert(y->size == dydt->size && y->size == 3);
    double S = gsl_vector_get(y, 0);
    double I = gsl_vector_get(y, 1);
    double R = gsl_vector_get(y, 2);
    double N = S+I+R;
    
    gsl_vector_set(dydt, 0, -I*S/N/Tc);
    gsl_vector_set(dydt, 1, I*S/N/Tc - I/Tr);
    gsl_vector_set(dydt, 2, I/Tr);
}
void (*sir_wrapper(double c, double r))(double, gsl_vector*,gsl_vector*){
      Tc = c; Tr = r;
      return func_sir;
}
int n = 3;
double *ms;
void func_nbody(double t, gsl_vector *y, gsl_vector *dydt){
    gsl_vector *rij = gsl_vector_alloc(3);
    gsl_vector_view ri;
    gsl_vector_view rj;
    gsl_vector_view ai;
    gsl_vector_view aj;
    //Format of y: x1, y1, z1, x2, y2, z2, ..., vx1, vy1, vz1, ...
    for(int i = 0; i<n; i++){
        for(int j = i+1; j<n; j++){
            //Calc acceleration
            ri = gsl_vector_subvector(y, 3*i, 3);
            rj = gsl_vector_subvector(y, 3*j, 3);
            gsl_vector_memcpy(rij, &rj.vector);
            gsl_blas_daxpy(-1, &ri.vector, rij);
            double d = gsl_blas_dnrm2(rij);
            
            ai = gsl_vector_subvector(dydt, 3*n + 3*i, 3);
            aj = gsl_vector_subvector(dydt, 3*n + 3*j, 3);
            
            double aij = ms[i]*ms[j]/(d*d*d);
            gsl_vector_scale(rij, aij);
            
            gsl_blas_daxpy(1, rij, &ai.vector);
            gsl_blas_daxpy(-1, rij, &aj.vector);
        }
    }
    gsl_vector_free(rij);
    
    //Update velocities
    gsl_vector_view v1 = gsl_vector_subvector(y, 3*n, 3*n);
    gsl_vector_view v2 = gsl_vector_subvector(dydt, 0, 3*n);
    gsl_vector_memcpy(&v2.vector, &v1.vector);
}
void (*nbody_wrapper(int _n, double *_ms))(double, gsl_vector*,gsl_vector*){
      n=_n; ms = _ms;
      return func_nbody;
}


int main(){
    // ------ SIMPLE HARMONIC MOTION ------ //
    printf("Simple harmonic motion as a test. u''=-u\nFrom 0 to 20 with max stepsize of 0.2\n");
    int maxsteps = 200, eqs = 2;
    gsl_matrix *ylist = gsl_matrix_alloc(maxsteps, eqs);
    gsl_vector *xlist = gsl_vector_alloc(maxsteps);
    
    //Initial value
    gsl_matrix_set(ylist, 0, 0, 1);
    gsl_matrix_set(ylist, 0, 1, 0);
    
    int steps = driver(func_shm, 0, 20, 0.001, 0.2, 1e-2, 1e-2, ylist, xlist);
    printf("Steps taken: %d\n\n", steps);
    
    FILE *SHM_out = fopen("SHM.out", "w");
    for(int i = 0; i<abs(steps); i++){
        gsl_vector_view yi = gsl_matrix_row(ylist, i);
        fprintf(SHM_out, "%g %g %g\n", gsl_vector_get(xlist, i), gsl_vector_get(&yi.vector, 0), gsl_vector_get(&yi.vector, 1));
    }
    gsl_matrix_free(ylist);
    gsl_vector_free(xlist);
    fclose(SHM_out);
    fprintf(stderr, "Simple harmonic motion done\n");
    
    // ------- SIR MODEL ------ //
    printf("\nSIR model\n");
    FILE *SIR_out = fopen("SIR.out", "w");
    maxsteps = 1000; eqs = 3;
    
    ylist = gsl_matrix_alloc(maxsteps, eqs);
    xlist = gsl_vector_alloc(maxsteps);
    
    //Initial value
    gsl_matrix_set(ylist, 0, 0, 1);
    gsl_matrix_set(ylist, 0, 1, 1e-2);
    gsl_matrix_set(ylist, 0, 2, 0);
    
    printf("(Tc, Tr) = (1, 7)\n");
    steps = driver(sir_wrapper(1, 7), 0, 100, 1e-7, 2e-1, 1e-4, 1e-4, ylist, xlist);
    printf("Steps taken: %d\n\n", steps);
    
    fprintf(SIR_out, "#index 0 - Tc=1, Tr=7\n");
    for(int i = 0; i<abs(steps); i++){
        gsl_vector_view yi = gsl_matrix_row(ylist, i);
        fprintf(SIR_out, "%g %g %g %g\n", gsl_vector_get(xlist, i), gsl_vector_get(&yi.vector, 0), gsl_vector_get(&yi.vector, 1), gsl_vector_get(&yi.vector, 2));
    }
    
    printf("(Tc, Tr) = (2, 7)\n");
    steps = driver(sir_wrapper(2, 7), 0, 100, 1e-7, 2e-1, 1e-4, 1e-4, ylist, xlist);
    printf("Steps taken: %d\n\n", steps);
    
    fprintf(SIR_out, "\n\n#index 1 - Tc=2, Tr=7\n");
    for(int i = 0; i<abs(steps); i++){
        gsl_vector_view yi = gsl_matrix_row(ylist, i);
        fprintf(SIR_out, "%g %g %g %g\n", gsl_vector_get(xlist, i), gsl_vector_get(&yi.vector, 0), gsl_vector_get(&yi.vector, 1), gsl_vector_get(&yi.vector, 2));
    }
    
    printf("(Tc, Tr) = (4, 7)\n");
    steps = driver(sir_wrapper(4, 7), 0, 100, 1e-7, 2e-1, 1e-4, 1e-4, ylist, xlist);
    printf("Steps taken: %d\n\n", steps);
    
    fprintf(SIR_out, "\n\n#index 2 - Tc=4, Tr=7\n");
    for(int i = 0; i<abs(steps); i++){
        gsl_vector_view yi = gsl_matrix_row(ylist, i);
        fprintf(SIR_out, "%g %g %g %g\n", gsl_vector_get(xlist, i), gsl_vector_get(&yi.vector, 0), gsl_vector_get(&yi.vector, 1), gsl_vector_get(&yi.vector, 2));
    }
    
    gsl_matrix_free(ylist);
    gsl_vector_free(xlist);
    
    fclose(SIR_out);
    fprintf(stderr, "Sir model done\n");
    
    // ------- NBODY ------ //
    printf("\nNBODY - N=3 and equal masses\n");
    FILE *nbody_out = fopen("nbody.out", "w");
    n = 3;
    double masses[3] = {1, 1, 1};
    maxsteps = 1000; eqs = 3*n*2;
    
    ylist = gsl_matrix_calloc(maxsteps, eqs);
    xlist = gsl_vector_alloc(maxsteps);
    
    //Initial value
    //Positions
    gsl_matrix_set(ylist, 0, 0, -0.97000436);
    gsl_matrix_set(ylist, 0, 1, 0.24308753);
    gsl_matrix_set(ylist, 0, 3, 0);
    gsl_matrix_set(ylist, 0, 4, 0);
    gsl_matrix_set(ylist, 0, 6, 0.97000436);
    gsl_matrix_set(ylist, 0, 7, -0.24308753);
    //Velocities
    gsl_matrix_set(ylist, 0, 3*n,   0.4662036850);
    gsl_matrix_set(ylist, 0, 3*n+1, 0.4323657300);
    gsl_matrix_set(ylist, 0, 3*n+3, -0.93240737);
    gsl_matrix_set(ylist, 0, 3*n+4, -0.86473146);
    gsl_matrix_set(ylist, 0, 3*n+6, 0.4662036850);
    gsl_matrix_set(ylist, 0, 3*n+7, 0.4323657300);
    
    steps = driver(nbody_wrapper(n, masses), 0, 5, 1e-3, 0.01, 1e-9, 1e-9, ylist, xlist);
    printf("Steps taken: %d\n\n", steps);
    
    for(int i = 0; i<abs(steps); i++){
        // fprintf(nbody_out, "%g ", gsl_vector_get(xlist, i));
        for(int j = 0; j<3*n; j++){
            fprintf(nbody_out, "%g ", gsl_matrix_get(ylist, i, j));
        }
        fprintf(nbody_out, "\n");
    }
    
    gsl_matrix_free(ylist);
    gsl_vector_free(xlist);
    
    fclose(nbody_out);
    return 0;
}