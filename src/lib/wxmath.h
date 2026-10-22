#ifndef __wxmath__h__
#define __wxmath__h__

#include <cmath>

#define SGN(b) ((b)>=0.)?1.:-1.

template <typename T>
T
dsqr(T a)
{ return a*a; }

template<typename T>
T
dmax(T a, T b)
{ return a>b?a:b; }

template<typename T>
T 
dmax(T a, T b, T c)
{ 
  T m1= a>b?a:b;
  return m1>c?m1:c;
}

template<typename T>
T 
dmin(T a, T b)
{ return a<b?a:b; }

template<typename T>
T 
dmin(T a, T b, T c)
{ 
  T m1= a<b?a:b;
  return m1<c?m1:c;
}

template<typename T>
T 
dsign(T a, T b)
{
  T sgn = (b>=0.0)?1.0:-1.0;
  return fabs(a)*sgn;
}

template<typename T>
T*
alloc_1d(unsigned nx, T val=0)
{
  T* r = new T[nx];
  for (unsigned i=0; i<nx; ++i)
    r[i] = val;
  return r;
}

template<typename T>
T**
alloc_2d_c(unsigned nx, unsigned ny, T val=0)
{
  // allocate a contigous piece of memory
  T *m = new T[nx*ny];
  // allocate an array of 'nx' pointers to pointers
  T **r = new T*[nx];
  // attach r to proper place in m
  for (unsigned i=0, j=0; i<nx; ++i, j=j+ny)
    r[i] = m+j;
  // initialize array
  for (unsigned i=0; i<nx; ++i)
    for (unsigned j=0; j<ny; ++j)
      r[i][j] = val;
  return r;    
}

template<typename T>
void
free_2d_c(T **r, unsigned nx, unsigned by)
{
  delete [] r[0]; // delete allocated memory
  delete [] r; // delete pointer to pointer
}

template<typename T>
T**
alloc_2d_p(unsigned nx, unsigned ny, T val=0)
{
  // allocate an array of 'nx' pointers to pointers
  T **r = new T*[nx];
  // allocate 'ny' Ts for each of the 'nx' pointers
  for(unsigned i=0; i<nx; ++i)
    r[i] = new T[ny];
  // initialize array
  for (unsigned i=0; i<nx; ++i)
    for (unsigned j=0; j<ny; ++j)
      r[i][j] = val;
  return r;    
}

template<typename T>
void
free_2d_p(T **r, unsigned nx, unsigned ny)
{
  // delete each of the 'nx' pointers
  for(unsigned i=0; i<nx; ++i)
    delete [] r[i];
  // delete the array of pointers to pointers
  delete [] r;
}

template<typename T>
T***
alloc_3d_c(unsigned nx, unsigned ny, unsigned nz, T val=0)
{
  T ***r = new T**[nx];
  r[0] = new T*[nx*ny];
  r[0][0] = new T[nx*ny*nz];
  for (unsigned j=1; j<ny; j++)
    r[0][j]=r[0][j-1] + nz;
  for (unsigned i=1; i<nx; i++)
  {
    r[i] = r[i-1] + ny;
    r[i][0] = r[i-1][0] + ny*nz;
    for (unsigned j=1; j<ny; j++)
      r[i][j] = r[i][j-1] + nz;
  }

  // initialize array
  for (unsigned i=0; i<nx; ++i)
    for (unsigned j=0; j<ny; ++j)
      for (unsigned k=0; k<nz; ++k)
        r[i][j][k] = val;
  return r; 
}

template<typename T>
void
free_3d_c(T ***r, unsigned nx, unsigned ny, unsigned nz)
{
  // delete each of the pointers
  if (r != 0) {
    delete[] (r[0][0]);
    delete[] (r[0]);
    delete[] (r);
  }
}

template<typename T>
void
swap_2d(T ***a, T ***b)
{
  T** t = *a;
  *a = *b;
  *b = t;
}

template<typename T>
void
swap_1d(T **a, T **b)
{
  T* t = *a;
  *a = *b;
  *b = t;
}

/**
   Compute absicca and weights of Gauss integration in the interval
   [x1,x2] of order n. The arrays x[1..n] and w[1..n] are set to the
   qudrature points and the weights repecively. The sum of the weights
   is normalized to 2.0.
*/
template<typename REAL>
void 
gauleg(REAL x1, REAL x2, REAL *x, REAL *w, int n)
{
  int m, j, i;
  REAL z1, z, xm, xl, pp, p3, p2, p1;
  REAL EPS = 3.0e-13;
  REAL PI = 3.14159265358979323846;

  m = (n+1)/2;
  xm = 0.5*(x2+x1);
  xl = 0.5*(x2-x1);
  for (i = 1; i <= m; i++)
  {
    z = cos(PI*(i-0.25)/(n+0.5));
    do
    {
      p1 = 1.0;
      p2 = 0.0;
      for (j = 1; j <= n; j++)
      {
        p3 = p2;
        p2 = p1;
        p1 = ((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
      }
      pp = n*(z*p1-p2)/(z*z-1.0);
      z1 = z;
      z = z1-p1/pp;
    } 
    while( fabs(z-z1) > EPS );
    x[i] = xm-xl*z;
    x[n+1-i] = xm+xl*z;
    w[i] = 2.0*xl/((1.0-z*z)*pp*pp);
    w[n+1-i] = w[i];
  }
}

/**
   Computes Legendre polynomial of order n at point x in [-1,1].
*/
template<class REAL>
REAL
legendre_p(int n, REAL x)
{
  REAL p0 = 1.0;
  REAL p1 = x;
  REAL pn,pn1,pn2;

  if (n==0) return p0;
  if (n==1) return p1;
  pn = 0.0; pn1 = p1; pn2 = p0; // initialize recurrence
  for(int i=2; i<=n; i++)
  {
    // use recurrence relation to compute P_n
    pn = (x*(2.*i-1.)*pn1 - (i-1.)*pn2)/(1.*i);
    pn2 = pn1;
    pn1 = pn;
  }
  return pn;
}

/**
   Computes derivative of the Legendre polynomial of order n at point
   x in [-1,1].
*/
template<class REAL>
REAL
legendre_p_d(int n, REAL x)
{
  REAL dp0 = 0;
  REAL dp1 = 1;
  REAL dpn;

  if (n==0) return dp0;
  if (n==1) return dp1;
  dpn = ((n+1.)*x*legendre_p(n,x) - (n+1.)*legendre_p(n+1,x))/(1.-x*x);

  return dpn;
}

/**
   Computes Lagrangian polynomial of order n at point x in [-1,1].
*/
template<class REAL>
REAL
lagrange_p(int n, REAL x, int num)
{
  REAL dx = 2.0/float(n);
  REAL p1 = 1.0;
  REAL XX[n+1];

  for(int i=0; i<n+1; ++i)
      XX[i]=-1.0+i*dx;

  for(int mn=0; mn<n+1; ++mn)
      if(mn!=num)
          p1 = p1*(x-XX[mn])/(XX[num]-XX[mn]);

  return p1;
}

/**
   Computes derivative of the Lagrange kth polynomial of order n at point
   x in [-1,1].
*/
template<class REAL>
REAL
lagrange_p_d(int n, REAL x, int k)
{
  REAL dx = 2.0/float(n);
  REAL p2 = 0.0;
  REAL XX[n+1];

  for(int i=0; i<n+1; ++i)
      XX[i]=-1.0+i*dx;

  for(int i=0; i<n+1; ++i)
  {
      if(i!=k)
      {
          REAL p1=1.0;
          for(int j=0;j<n+1; ++j)
              if((i!=j) && (j!=k))
                  p1 = p1*((x-XX[j])/(XX[k]-XX[j]));

          p2 += p1/(XX[k]-XX[i]);
      }
  }

  return p2;
}

/** Jacobi polynomials */
template<class REAL>
REAL
JacobiP(REAL x, REAL alpha, REAL beta, int N)
//---------------------------------------------------------
{
  // function [P] = JacobiP(x,alpha,beta,N)
  // Purpose: Evaluate Jacobi Polynomial of type (alpha,beta) > -1
  //          (alpha+beta <> -1) at point x for order N and
  //          returns P
  // Note   : They are normalized to be orthonormal.

  REAL aold=0.0, anew=0.0, bnew=0.0, h1=0.0;
  REAL gamma0=0.0, gamma1=0.0;
  REAL ab=alpha+beta, ab1=alpha+beta+1.0, a1=alpha+1.0, b1=beta+1.0;

  // PL2 is assigned by the recurrence below, which runs only for N >= 2;
  // initialise it so an out-of-range N returns a value rather than reading
  // an indeterminate one.
  REAL P = 0.0, PL0 = 0.0, PL1 = 0.0, PL2 = 0.0, prow = 0.0, x_bnew = 0.0;

  // Initial values P_0(x) and P_1(x)
  gamma0 = pow(2.0,ab1)/(ab1)*tgamma(a1)*tgamma(b1)/tgamma(ab1);

  if (0==N) { P   = 1.0/sqrt(gamma0);  return P;
  } else { PL0 = 1.0/sqrt(gamma0); }

  gamma1 = (a1)*(b1)/(ab+3.0)*gamma0;
  prow = ((ab+2.0)*x/2.0 + (alpha-beta)/2.0) / sqrt(gamma1);

  if (1==N) { P   = prow; return P;
  } else { PL1 = prow; }

  // Repeat value in recurrence.
  aold = 2.0/(2.0+ab)*sqrt((a1)*(b1)/(ab+3.0));

  // Forward recurrence using the symmetry of the recurrence.
  for (int i=1; i<=(N-1); ++i) {
    h1 = 2.0*i+ab;
    anew = 2.0/(h1+2.0)*sqrt((i+1)*(i+ab1)*(i+a1)*(i+b1)/(h1+1.0)/(h1+3.0));
    bnew = - (alpha*alpha-beta*beta)/h1/(h1+2.0);
    x_bnew = x-bnew;
    PL2 = 1.0/anew*( -aold*PL0 + x_bnew*PL1);
    aold =anew; PL0 = PL1; PL1 = PL2;
  }

  P = PL2;
  return P;
}

/** Jacobi polynomials */
template<class REAL>
REAL
GradJacobiP(REAL z,REAL alpha,REAL beta,int N)
//---------------------------------------------------------
{
  // function [dP] = GradJacobiP(z, alpha, beta, N);
  // Purpose: Evaluate the derivative of the orthonormal Jacobi
  //	   polynomial of type (alpha,beta)>-1, at point x
  //          for order N and returns dP

  REAL dP;
  if (0 == N) {
    dP=0.0;
  } else {
    dP = sqrt(N*(N+alpha+beta+1))*JacobiP(z,alpha+1,beta+1, N-1);
  }
  return dP;
}


template<class REAL>
REAL
mminmod(REAL a, REAL b, REAL c, REAL dx, REAL M)
{
  if (fabs(a) < M*dx*dx)
    return a;
  double sa = SGN(a);
  double sb = SGN(b);
  double sc = SGN(c);
  if( (sa==sb) && (sb==sc) )
  {
    if (sa<0)
      return dmax(a,b,c);
    else
      return dmin(a,b,c);
  }
  else
    return 0;
}

template<class REAL>
void
frobsum(REAL &frobsum, REAL mat1[3][3], REAL mat2[3][3])
{
	for (int m=0; m<3; ++m)
		for (int nn=0; nn<3; ++nn)
			frobsum      += mat1[m][nn]*mat2[m][nn];

//	return frobsum;
}

template<class REAL>
void
vectorcrossmatrix(REAL vecxmat[3][3], REAL vec[3], REAL mat[3][3])
{
	int vld, mld, vlg, mlg;

	for (int jj=0; jj<3; ++jj)
		for (int ii=0; ii<3; ++ii)
		{
			if (ii==0)
			{
				vld = 1;
				mld = 2;
				vlg = 2;
				mlg = 1;
			}
			else if (ii==1)
			{
				vld = 2;
				mld = 0;
				vlg = 0;
				mlg = 2;
			}
			else if (ii==2)
			{
				vld = 0;
				mld = 1;
				vlg = 1;
				mlg = 0;
			}
			vecxmat[ii][jj]= vec[vld]*mat[mld][jj]-vec[vlg]*mat[mlg][jj];
		}

//	return vecxmat;
}


template<class REAL>
void
vectorcrossvector(REAL vecxvec[3], REAL vec1[3], REAL vec2[3])
{
	int vld, mld, vlg, mlg;

	for (int ii=0; ii<3; ++ii)
	{
		if (ii==0)
		{
			vld = 1;
			mld = 2;
			vlg = 2;
			mlg = 1;
		}
		else if (ii==1)
		{
			vld = 2;
			mld = 0;
			vlg = 0;
			mlg = 2;
		}
		else if (ii==2)
		{
			vld = 0;
			mld = 1;
			vlg = 1;
			mlg = 0;
		}
		vecxvec[ii]= vec1[vld]*vec2[mld]-vec1[vlg]*vec2[mlg];
	}
//	return vecxvec;
}


template<class REAL>
void
vectordotvector(REAL vecdotvec[3], REAL vec1[3], REAL vec2[3])
{
	vecdotvec = vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2];
//	return vecdotvec;
}



template<class REAL>
void
vecvec(REAL mat[3][3], REAL vec1[3], REAL vec2[3])
{
	for (int i=0; i<3; ++i)
		for (int j=0; j<3; ++j)
			mat[i][j] = vec1[i]*vec2[j];
//	returns tensor;
}


template<class REAL>
void
vectordotmatrix(REAL vecdotmat[3], REAL vec[3], REAL mat[3][3])
{
	for (int i=0; i<3; ++i)
		vecdotmat[i] = vec[0]*mat[0][i] + vec[1]*mat[1][i] + vec[2]*mat[2][i];
//	return vecdotmat;
}

template<class REAL>
void
normalizevector(REAL vec[3])
{
    REAL vecLen = sqrt(vec[0]*vec[0]
      + vec[1]*vec[1] + vec[2]*vec[2]);

    // make sure it is unit tangent
    vec[0] = vec[0]/vecLen;
    vec[1] = vec[1]/vecLen;
    vec[2] = vec[2]/vecLen;

//    return vec;
}

#endif // __wxmath__h__
