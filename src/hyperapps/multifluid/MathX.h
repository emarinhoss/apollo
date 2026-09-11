/*
 * File:   MathX.h
 * Author: cambier
 *
 * Created on July 21, 2011, 2:00 PM
 */

#ifndef MATHX_H
#define	MATHX_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cmath>
#include <time.h>
#include <vector>
#include "constant.h"
#include "kmatrix.h"
#include <Eigen/Dense>  // Eigen3 library (eigen3 dir must be in include path)

using namespace std;

class MathX
{
private:
public:
    MathX() { };
   ~MathX() { };


    static double gamma( double );    // Gamma function
    static double erf ( double );     // error function
    static double erfc( double );     // complement to error function
    static double En(int, double , double);  // exponential integral (int_a^b(dy/y^n e^-y)
    static double Psi(double);
    static double expIntegral(int,double);

    // Legendre polynomials
    static double** LegendrePolyCoefs(int order);
    static double  LegendrePolynomial(int, double*, double);
    static double* LegendrePolynomial(int, double*, int, double*);

    // Gaussian quadrature
    static kMatrix2D<double>* gaussQuadNodes(int);
    static kMatrix2D<double>* gaussQuadNodes(int, double, double);

    // splines
    static double* getSplineCoefs(int, double*, double* );
    static double* computeSpline(int, double*, double*, double*);

    // interpolation
    template<class T> static T interpolate1d(T* xb,T x,T* dat)
    {
    	T xd = (x-xb[0])/(xb[1]-xb[0]);
    	T c   = dat[0]  * (1-xd) + dat[1]  * xd;
    	return c;
    }

    template<class T> static T interpolate2d(T* xb,T* yb,T x,T y,T** dat)
    {
    	T xd = (x-xb[0])/(xb[1]-xb[0]);
    	T yd = (y-yb[0])/(yb[1]-yb[0]);
    	T c0  = dat[0][0] * (1-xd) + dat[1][0] * xd;
    	T c1  = dat[0][1] * (1-xd) + dat[1][1] * xd;
    	T c   = c0  * (1-yd) + c1  * yd;
    	return c;
    }

    template<class T> static T interpolate3d(T* xb,T* yb,T* zb,T x,T y,T z,T*** dat)
    {
    	T xd = (x-xb[0])/(xb[1]-xb[0]);
    	T yd = (y-yb[0])/(yb[1]-yb[0]);
    	T zd = (z-zb[0])/(zb[1]-zb[0]);
    	T c00 = dat[0][0][0] * (1.-xd) + dat[1][0][0] * xd;
    	T c10 = dat[0][1][0] * (1.-xd) + dat[1][1][0] * xd;
    	T c01 = dat[0][0][1] * (1.-xd) + dat[1][0][1] * xd;
    	T c11 = dat[0][1][1] * (1.-xd) + dat[1][1][1] * xd;
    	T c0  = c00 * (1-yd) + c10 * yd;
    	T c1  = c01 * (1-yd) + c11 * yd;
    	T c   = c0  * (1-zd) + c1  * zd;
    	return c;
    }


/**
 *  Gauss-Jordan elimination (implicit point-solver)
 *  Solves A*s = r for s
 * @param int N = size
 * @param s[]  = solution (vector of size N)
 * @param A[][]= matrix size NxN
 * @param r[]  = right-hand-side (vector of size N)
 */
    static void Gauss_Jordan_point(int, float* , float** , float* );
    static void Gauss_Jordan_point(int, double*, double**, double*);
    static void EigenDecomp(int, double**, double**, double**, double*);
    static void EigenDecomp(int, double**, complex<double>**, complex<double>**, complex<double>*);
    static void LU_decomp (int, double**, int*);
    static double** invert(int, double**);


};

class LegendrePolynomials
{
private:
    int Np;
    double** cp;
public:
    LegendrePolynomials() { }

    void setCoefficients(int maxOrder);
    double evalPolynomial(int order, double z);

}; // end class LegendrePolynomials

class GaussQuadrature
{
private:
public:
    int Nq;
    double* qcoo;
    double* qwgt;

    GaussQuadrature(int);
    GaussQuadrature(int,double,double);
   ~GaussQuadrature();

    void quadNodes(int,double,double);

};

//// C-array utilities
//using namespace std;
//// 1D
//template<class T>  T* cArray(int N0)
//{
//    T* A = new T[N0];
//    for(int i=0;i<N0;i++) A[i] = (T)0;
//    return A;
//}
//template<class T> void del_cArray(int N0, T* A)
//{
//    if(A)
//    {
//        delete[] A;
//    }
//}
//
//// 2D
//template<class T>  T** cArray(int N0, int N1)
//{
//    T** A = new T*[N0];
//    for(int i=0;i<N0;i++)
//    {
//        A[i] = new T[N1];
//        for(int j=0;j<N1;j++) A[i][j] = (T)0;
//    }
//    return A;
//}
//template<class T> void del_cArray(int N0, int N1, T** A)
//{
//    if(A)
//    {
//        for(int i=0;i<N0;i++) delete[] A[i];
//        delete[] A;
//    }
//}
//
//// 3D
//template<class T>  T*** cArray(int N0, int N1, int N2)
//{
//    T*** A = new T**[N0];
//    for(int i=0;i<N0;i++)
//    {
//        A[i] = new T*[N1];
//        for(int j=0;j<N1;j++)
//        {
//            A[i][j] = new T[N2];
//            for(int k=0;k<N2;k++) A[i][j][k] = (T)0;
//        }
//    }
//    return A;
//}
//template<class T> void del_cArray(int N0, int N1, int N2, T*** A)
//{
//    if(A)
//    {
//        for(int i=0;i<N0;i++)
//        {
//            for(int j=0;j<N1;j++) delete[] A[i][j];
//            delete[] A[i];
//        }
//        delete[] A;
//    }
//}
//
//// 4D
//template<class T>  T**** cArray(int N0, int N1, int N2, int N3)
//{
//    T**** A = new T***[N0];
//    for(int i=0;i<N0;i++)
//    {
//        A[i] = new T**[N1];
//        for(int j=0;j<N1;j++)
//        {
//            A[i][j] = new T*[N2];
//            for(int k=0;k<N2;k++)
//            {
//                A[i][j][k] = new T[N3];
//                for(int m=0;m<N3;m++) A[i][j][k][m] = (T)0;
//            }
//        }
//    }
//    return A;
//}
//template<class T> void del_cArray(int N0, int N1, int N2, int N3, T**** A)
//{
//    if(A)
//    {
//        for(int i=0;i<N0;i++)
//        {
//            for(int j=0;j<N1;j++)
//            {
//                for(int k=0;k<N2;k++) delete[] A[i][j][k];
//                delete[] A[i][j];
//            }
//            delete[] A[i];
//        }
//        delete[] A;
//    }
//}
//
//// 5D
//template<class T>  T***** cArray(int N0, int N1, int N2, int N3, int N4)
//{
//    T***** A = new T****[N0];
//    for(int i=0;i<N0;i++)
//    {
//        A[i] = new T***[N1];
//        for(int j=0;j<N1;j++)
//        {
//            A[i][j] = new T**[N2];
//            for(int k=0;k<N2;k++)
//            {
//                A[i][j][k] = new T*[N3];
//                for(int m=0;m<N3;m++)
//                {
//                    A[i][j][k][m] = new T[N4];
//                    for(int n=0;n<N4;n++) A[i][j][k][m][n] = (T)0;
//                }
//            }
//        }
//    }
//    return A;
//}
//template<class T> void del_cArray(int N0, int N1, int N2, int N3, int N4, T***** A)
//{
//    if(A)
//    {
//        for(int i=0;i<N0;i++)
//        {
//            for(int j=0;j<N1;j++)
//            {
//                for(int k=0;k<N2;k++)
//                {
//                    for(int m=0;m<N3;m++) delete[] A[i][j][k][m];
//                    delete[] A[i][j][k];
//                }
//                delete[] A[i][j];
//            }
//            delete[] A[i];
//        }
//        delete[] A;
//    }
//}
//
//// 6D
//template<class T>  T****** cArray(int N0, int N1, int N2, int N3, int N4, int N5)
//{
//    T****** A = new T*****[N0];
//    for(int i=0;i<N0;i++)
//    {
//        A[i] = new T****[N1];
//        for(int j=0;j<N1;j++)
//        {
//            A[i][j] = new T***[N2];
//            for(int k=0;k<N2;k++)
//            {
//                A[i][j][k] = new T**[N3];
//                for(int m=0;m<N3;m++)
//                {
//                    A[i][j][k][m] = new T*[N4];
//                    for(int n=0;n<N4;n++)
//                    {
//                        A[i][j][k][m][n] = new T[N5];
//                        for(int p=0;p<N5;p++) A[i][j][k][m][n][p] = (T)0;
//                    }
//                }
//            }
//        }
//    }
//    return A;
//}
//template<class T> void del_cArray(int N0, int N1, int N2, int N3, int N4, int N5, T****** A)
//{
//    if(A)
//    {
//        for(int i=0;i<N0;i++)
//        {
//            for(int j=0;j<N1;j++)
//            {
//                for(int k=0;k<N2;k++)
//                {
//                    for(int m=0;m<N3;m++)
//                    {
//                        for(int n=0;n<N4;n++) delete[] A[i][j][k][m][n];
//                        delete[] A[i][j][k][m];
//                    }
//                    delete[] A[i][j][k];
//                }
//                delete[] A[i][j];
//            }
//            delete[] A[i];
//        }
//        delete[] A;
//    }
//}
//
//// 7D
//template<class T>  T******* cArray(int N0, int N1, int N2, int N3, int N4, int N5, int N6)
//{
//    T******* A = new T******[N0];
//    for(int i=0;i<N0;i++)
//    {
//        A[i] = new T*****[N1];
//        for(int j=0;j<N1;j++)
//        {
//            A[i][j] = new T****[N2];
//            for(int k=0;k<N2;k++)
//            {
//                A[i][j][k] = new T***[N3];
//                for(int m=0;m<N3;m++)
//                {
//                    A[i][j][k][m] = new T**[N4];
//                    for(int n=0;n<N4;n++)
//                    {
//                        A[i][j][k][m][n] = new T*[N5];
//                        for(int p=0;p<N5;p++)
//                        {
//                            A[i][j][k][m][n][p] = new T[N6];
//                            for(int q=0;q<N6;q++) A[i][j][k][m][n][p][q] = (T)0;
//                        }
//                    }
//                }
//            }
//        }
//    }
//    return A;
//}
//template<class T> void del_cArray(int N0, int N1, int N2, int N3, int N4, int N5, int N6, T******* A)
//{
//    if(A)
//    {
//        for(int i=0;i<N0;i++)
//        {
//            for(int j=0;j<N1;j++)
//            {
//                for(int k=0;k<N2;k++)
//                {
//                    for(int m=0;m<N3;m++)
//                    {
//                        for(int n=0;n<N4;n++)
//                        {
//                            for(int p=0;p<N5;p++) delete[] A[i][j][k][m][n][p];
//                            delete[] A[i][j][k][m][n];
//                        }
//                        delete[] A[i][j][k][m];
//                    }
//                    delete[] A[i][j][k];
//                }
//                delete[] A[i][j];
//            }
//            delete[] A[i];
//        }
//        delete[] A;
//    }
//}


#endif	/* MATHX_H */



