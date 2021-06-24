#include"matrix.h"
#include<gsl/gsl_blas.h>
#include<math.h>
#include<assert.h>

//Method for calculating matrix product between some square matrix A and a rotation matrix in O(n)
// A <- AJ(theta)
void timesJ(gsl_matrix *A, int p, int q, double theta){
    double c = cos(theta), s = sin(theta);
    double Aip, Aiq;
    //Changes the columns p and q in A
    for(int i = 0; i<A->size1; i++){
        Aip = gsl_matrix_get(A, i, p);
        Aiq = gsl_matrix_get(A, i, q);
        gsl_matrix_set(A, i, p, c*Aip-s*Aiq);
        gsl_matrix_set(A, i, q, s*Aip+c*Aiq);
    }
}
//Method for calculating matrix product between a rotation matrix and some square matrix A in O(n)
// A <- J(theta)A
void Jtimes(gsl_matrix *A, int p, int q, double theta){
    double c = cos(theta), s = sin(theta);
    double Apj, Aqj;
    //Changes the rows p and q in A
    for(int j = 0; j<A->size2; j++){
        Apj = gsl_matrix_get(A, p, j);
        Aqj = gsl_matrix_get(A, q, j);
        gsl_matrix_set(A, p, j,  c*Apj + s*Aqj);
        gsl_matrix_set(A, q, j, -s*Apj + c*Aqj);
    }
}

int SVD_two_jaco_square(gsl_matrix *A, gsl_matrix *V, gsl_matrix *U){
    assert(A->size1==A->size2);
    int n = A->size1;

    int changes = 0;
    int sweeps = 0;

    double old_app, old_aqq, new_app, new_aqq, apq, aqp, app, aqq;
    double givens_angle, jacobi_angle;
    //Reset U,V to identity
    gsl_matrix_set_identity(V);
    gsl_matrix_set_identity(U);
    do{
        changes=0;
        sweeps++;
        //Sweep over upper triangle
        for(int p=0;p<n-1;p++){
            for(int q=p+1;q<n;q++){
                //Check if diagonal elements change upon next rotation
                old_app=gsl_matrix_get(A,p,p);
                old_aqq=gsl_matrix_get(A,q,q);

                //Calculate A <- J^T G^T A J 
                apq=gsl_matrix_get(A,p,q);
                aqp=gsl_matrix_get(A,q,p);
                app=gsl_matrix_get(A,p,p);
                aqq=gsl_matrix_get(A,q,q);

                //First symmetrizise the two off diagonal elements
                givens_angle=atan2(apq-aqp, app+aqq);            //Angle for Givens matrix
                Jtimes(A,p,q, -givens_angle);

                apq=gsl_matrix_get(A,p,q);
                app=gsl_matrix_get(A,p,p);
                aqq=gsl_matrix_get(A,q,q);

                //Then eliminate them
                jacobi_angle=0.5*atan2(2*apq,aqq-app);            //Angle for Jacobi matrix
                timesJ(A,p,q,  jacobi_angle);
                Jtimes(A,p,q, -jacobi_angle);

                //Update U and V
                timesJ(V,p,q, jacobi_angle);
                timesJ(U,p,q, givens_angle);
                timesJ(U,p,q, jacobi_angle);

                new_app=gsl_matrix_get(A,p,p);
                new_aqq=gsl_matrix_get(A,q,q);

                //Convergence criteria - the new diagonal elements must be equal to the old, after one sweep (within machine epsilon)
                if(old_app != new_app || old_aqq != new_aqq){
                    changes++;
                }
            }
        }
    } while(changes>0);

    //Loop through and check for negative diagonal elements in D
    for(int i = 0; i<n; i++){
        gsl_vector_view Ucol;
        double dii = gsl_matrix_get(A, i, i);
        if(dii < 0){
            Ucol = gsl_matrix_column(U, i);
            gsl_vector_scale(&Ucol.vector, -1);
            gsl_matrix_set(A, i, i, -dii);
        }
    }
    return sweeps;
}

//SVD on general tall matrix - A is destroyed in the process.
//Or rather transformed in to Q of the QR factorization of A.
int SVD_two_jaco(gsl_matrix *A, gsl_matrix *D, gsl_matrix *V, gsl_matrix *U){
    assert(A->size1>A->size2);
    
    //First QR-factorize A = Q*R, with R=D
    GS_decomp(A, D);
    
    //SVD on square matrix R=U*R*VT
    gsl_matrix *Uprime = gsl_matrix_alloc(A->size2, A->size2);
    int sweeps = SVD_two_jaco_square(D, V, Uprime);
    
    //Transforming U <- Q*U'
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, A, Uprime, 0, U);
    
    gsl_matrix_free(Uprime);
    return sweeps;
}