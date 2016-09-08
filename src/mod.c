#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#define MAT2(mat, i, j, nI) (mat[i+j*nI])
#define MAT3(mat, i, j, a, nI, nJ) (mat[i+j*nI+a*nI*nJ])



void modification(
                    int *p,
                    double *InfoInv,
                    double *nu_stu,
                    double *nu_s_tu,
                    double *mod
                    )
{
    double cons;
    int r,s,t,u,P=p[0];
    
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