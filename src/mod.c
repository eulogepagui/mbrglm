#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
//#include <R_ext/Linpack.
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#define MAT2(mat, i, j, nI) (mat[i+j*nI])
#define MAT3(mat, i, j, a, nI, nJ) (mat[i+j*nI+a*nI*nJ])


void Matxvet(int nrow,  //dim is the length of vet
             int ncol,
             double *vet,
             double *matIn,
             double *matOut
             )
{
    int i,j;
    for(i=0; i <nrow; i++)
    {
        for(j=0; j <ncol; j++)
        {
            matOut[i + nrow*j]= matIn[i + nrow*j]*vet[i];
        }
    }
}


//--- product of 2 matrices A^T%*%B
void prod2mat2(
               int nrowA,
               int ncolB,
               int ncolA,
               double alpha,
               double *a,
               double *b,
               double *ab
               )
{  
    char transA='t';
    char transB='n';
    double beta=0.0; 
    F77_CALL(dgemm)(&transA, &transB, &ncolA,&ncolB, &nrowA, &alpha,a, &nrowA,b, &nrowA,&beta, ab, &ncolA);
}  



void modification(
                  int *p,
                  int *n,
                  double *X,
                  double *mueta,
                  double *mu,
                  double *A,
                  double *B,
                  double *InfoInv,
                  double *weights,
                  double *mod
                    )
{
    double alpha=1.0,cons;
    int i,r,s,t,u,P=p[0],N=n[0],NN=N*N,NP=N*P,PP=P*P,PPP=P*P*P;
    double *nu_stu = (double *)  R_alloc(PPP, sizeof(double));
    double *nu_s_tu = (double *)  R_alloc(PPP, sizeof(double));
    double *vet = (double *)  R_alloc(NN, sizeof(double));
    double *vet2 = (double *)  R_alloc(NN, sizeof(double));
    double *mat = (double *)  R_alloc(NP, sizeof(double));
    double *matarray;
    double *matarray2;
    matarray=nu_stu;
    matarray2=nu_s_tu;
    
    for (r=0; r<P; r++)
    {
        for (i=0; i<N; i++) {
            vet[i]=R_pow_di(A[i],3)*mu[i]*(1-mu[i])*(1-2*mu[i])*MAT2(X,i,r,N)/R_pow_di(weights[i],2);
            vet2[i] = B[i]*mueta[i]*MAT2(X,i,r,N);
        }
        Matxvet(N,P,vet,X,mat);
        prod2mat2(N,P,P,alpha,X,mat,matarray);
        matarray=nu_stu+r*PP;
        matarray+=PP;
        Matxvet(N,P,vet2,X,mat);
        prod2mat2(N,P,P,alpha,X,mat,matarray2);
        matarray2=nu_s_tu+r*PP;
        matarray2+=PP;
    }

    
    for (r=0; r<P; r++)
    {
        for (s=0; s<P; s++)
        {
            for (t=0; t<P; t++)
            {
                for (u=t; u<P; u++)
                {
                    if(t==u)cons=1.0;
                    else cons=2.0;
                    mod[r] = (mod[r] +cons*0.5*MAT2(InfoInv,r,s,P)*(MAT2(InfoInv,t,u,P)*(MAT3(nu_stu,s,t,u,P,P)+MAT3(nu_s_tu,s,t,u,P,P))
							-MAT2(InfoInv,r,t,P)*MAT2(InfoInv,r,u,P)*(2.0*MAT3(nu_stu,s,t,u,P,P)/3.0+MAT3(nu_s_tu,s,t,u,P,P))/MAT2(InfoInv,r,r,P)) );
                }
            }
        }

    }
}	