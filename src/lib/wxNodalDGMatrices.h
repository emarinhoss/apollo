#ifndef WXNODALDGMATRICES_H
#define WXNODALDGMATRICES_H

#include <petsc.h>
#include "wxmath.h"

template<typename REAL>
REAL
Simplex2DP(REAL a, REAL b, int i, int j)
{
    // function [P] = Simplex2DP(a,b,i,j);
    // Purpose : Evaluate 2D orthonormal polynomial
    //           on simplex at (a,b) of order (i,j).

    REAL P, number=2.*i+1., zero=0;
    REAL h1 = JacobiP(a,zero,zero,i), h2 = JacobiP(b,number,zero,j);
    REAL tv1=sqrt(2.)*h1*h2, tv2=pow(1.0-b,(REAL)i);
    P = tv1*tv2;
    return P;
}

//---------------------------------------------------------
template<typename REAL>
void rstoab(REAL r, REAL s, REAL aa, REAL bb)
//---------------------------------------------------------
{
    // function [a,b] = rstoab(r,s)
    // Purpose : Transfer from (r,s) -> (a,b) coordinates in triangle

    bb = -1.;
    if (s == 1.0){aa = -1.;}
    else{aa = 2.0*(1.+r)/(1.0-s)-1.0;}
}

//---------------------------------------------------------
template<typename REAL>
void
Vandermonde2D(int N, int Nx, REAL *rr, REAL *ss, Mat *Vand)
//---------------------------------------------------------
{
  // function [V2D] = Vandermonde2D(N, r, s);
  // Purpose : Initialize the 2D Vandermonde Matrix.
  //           V_{ij} = phi_j(r_i, s_i);

  REAL a,b;
  // build the Vandermonde matrix
  int sk = 0;
  for (int i=0; i<=N; ++i) {
    for (int j=0; j<=(N-i); ++j) {
        for (int k=0; k<Nx;k++)
        {
            REAL r = rr[k], s = ss[k];
            b = s;
            if (s == 1.0){a = -1.;}
            else{a = 2.0*(1.+r)/(1.0-s)-1.0;}
            PetscScalar value = Simplex2DP(a,b,i,j);
            MatSetValue(*Vand,k,sk,value,INSERT_VALUES);
        }
        ++sk;
    }
  }

  MatAssemblyBegin(*Vand,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*Vand,MAT_FINAL_ASSEMBLY);
}

//---------------------------------------------------------
template<typename REAL>
void GradSimplex2DP(const REAL a, const REAL b, int id, int jd, REAL *pdmodedr, REAL *pdmodeds)
//---------------------------------------------------------
{
  // function [dmodedr, dmodeds] = GradSimplex2DP(a,b,id,jd)
  // Purpose: Return the derivatives of the modal basis (id,jd)
  //          on the 2D simplex at (a,b).

  REAL fa, dfa, gb, dgb, tmp;
  REAL &dmodedr = *pdmodedr;
  REAL &dmodeds = *pdmodeds;

  fa = JacobiP(a, (REAL)0.0, (REAL)0.0, id);     dfa = GradJacobiP(a, (REAL)0.0, (REAL)0.0, id);
  gb = JacobiP(b, (REAL)2.0*id+1,(REAL)0.0, jd); dgb = GradJacobiP(b,(REAL)2.0*id+1,(REAL)0.0, jd);

  // r-derivative
  // d/dr = da/dr d/da + db/dr d/db = (2/(1-s)) d/da = (2/(1-b)) d/da
  dmodedr = dfa*gb;
  if (id>0) {
    dmodedr *= pow(0.5*(1.0-b), (id-1.0));
  }

  // s-derivative
  // d/ds = ((1+a)/2)/((1-b)/2) d/da + d/db
  dmodeds = dfa*gb*0.5*(1.0+a);
  if (id>0) {
    dmodeds *= pow(0.5*(1.0-b), (id-1.0));
  }

  tmp = dgb*pow(0.5*(1.0-b), REAL(id));
  if (id>0) {
    tmp -= (0.5*id)*gb*pow(0.5*(1.0-b), (id-1.0));
  }
  dmodeds += fa*tmp;

  // Normalize
  dmodedr *= pow(2.0, (id+0.5)); dmodeds *= pow(2.0, (id+0.5));
}

//---------------------------------------------------------
template<typename REAL>
void GradVandermonde2D(int   N, int Nr, REAL *rr, REAL *ss,  Mat V2Dr, Mat V2Ds)
//---------------------------------------------------------
{
  // function [V2Dr,V2Ds] = GradVandermonde2D(N,r,s)
  // Purpose : Initialize the gradient of the modal basis (i,j)
  //		at (r,s) at order N

    // Transfer to (a,b) coordinates
    REAL a=0., b=0., ddr=0., dds=0.;

    // build the Vandermonde matrix
    int sk = 0;
    for (int i=0; i<=N; ++i)
      for (int j=0; j<=(N-i); ++j) {
          for (int k=0; k<Nr;k++)
          {
              REAL r = rr[k], s = ss[k];
              b = s;
              if (s == 1.0){a = -1.;}
              else{a = 2.0*(1.+r)/(1.0-s)-1.0;}
              GradSimplex2DP(a,b,i,j, &ddr,&dds);
              MatSetValue(V2Dr,k,sk,ddr,INSERT_VALUES);
              MatSetValue(V2Ds,k,sk,dds,INSERT_VALUES);
          }
          ++sk;
      }

    MatAssemblyBegin(V2Dr,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(V2Dr,MAT_FINAL_ASSEMBLY);

    MatAssemblyBegin(V2Ds,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(V2Ds,MAT_FINAL_ASSEMBLY);

}

//---------------------------------------------------------
template<typename REAL>
void DifferentiationMatrices2D
(int   N, int Nr, REAL *r, REAL *s, Mat IVand, Mat *Dr, Mat *Ds)
//---------------------------------------------------------
{
  // function [Dr,Ds] = Dmatrices2D(N,r,s,V)
  // Purpose : Initialize the (r,s) differentiation matrices
  //	    on the simplex, evaluated at (r,s) at order N

  int Nodes = (N+1)*(N+2)/2;
  Mat Vr, Vs;
  MatCreateSeqDense(PETSC_COMM_SELF,Nr,Nodes,PETSC_NULL,&Vr);
  MatCreateSeqDense(PETSC_COMM_SELF,Nr,Nodes,PETSC_NULL,&Vs);


  GradVandermonde2D(N, Nr, r, s, Vr, Vs);
  //MatView(Vr,PETSC_VIEWER_STDOUT_WORLD);
  //MatView(Vs,PETSC_VIEWER_STDOUT_WORLD);
  MatMatMult(Vr,IVand,MAT_INITIAL_MATRIX,PETSC_DEFAULT,Dr);
  MatMatMult(Vs,IVand,MAT_INITIAL_MATRIX,PETSC_DEFAULT,Ds);

  MatDestroy(&Vr);
  MatDestroy(&Vs);
}

#endif // WXNODALDGMATRICES_H
