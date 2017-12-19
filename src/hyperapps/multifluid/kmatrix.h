/*
 * File:   kmatrix.h
 * Author: cambier
 *
 * Created on December 6, 2011, 2:01 PM
 * Revised on  August 27, 2011, 3:37 PM
 * Revised on January 15, 2012, 9:23 AM
 */

#ifndef KMATRIX_H
#define	KMATRIX_H

#ifndef NVCC
    #define NVCC 0
#endif

#if NVCC
    #define HDV __host__ __device__
    #define HST __host__
    #define DEV __device__
#else
    #define HDV
    #define HST
    #define DEV
    typedef struct { float x; float y; float z; float w; } float4;
#endif

#include "constant.h"
#include <stdio.h>
#include <stdlib.h>

#ifndef MATH_OPTS_KMAT
    #define MATH_OPTS_KMAT
    #if FOR_ARRAY
        #define kMatrix1D kFArray1D
        #define kMatrix2D kFArray2D
        #define kMatrix3D kFArray3D
        #define kMatrix4D kFArray4D
        #define kMatrix5D kFArray5D
        #define kMatrix6D kFArray6D
        #define kMatrix7D kFArray7D
        #define kF4Matrix1D kF4Array1D
    #else
        #define kMatrix1D kCArray1D
        #define kMatrix2D kCArray2D
        #define kMatrix3D kCArray3D
        #define kMatrix4D kCArray4D
        #define kMatrix5D kCArray5D
        #define kMatrix6D kCArray6D
        #define kMatrix7D kCArray7D
    #endif
#endif

typedef struct { int i0; }                   index1D;
typedef struct { int i0,i1; }                index2D;
typedef struct { int i0,i1,i2; }             index3D;
typedef struct { int i0,i1,i2,i3; }          index4D;
typedef struct { int i0,i1,i2,i3,i4; }       index5D;
typedef struct { int i0,i1,i2,i3,i4,i5; }    index6D;
typedef struct { int i0,i1,i2,i3,i4,i5,i6; } index7D;

/**
 * Matrices stored in Row-Major format have last dimension varying first:
 * for example, the indices of the elements of a 3D array of
 * dimension IMAX x JMAX x KMAX in row-major order are given by
 *       i*JMAX*KMAX + j*KMAX + k
 *
 * Matrices stored in Column-Major format have first dimension varying first:
 * for example, the indices of the elements of a 3D array of
 * dimension IMAX x JMAX x KMAX in col-major order are given by
 *       k*JMAX*IMAX+j*IMAX + i
 *
 * The default for Fortran arrays is Column-Major Format
 * The default for C-arrays (arrays of arrays) is Row-major format
 */


template<class type> class kFArray1D
{ // col-major storage: first index varies first
public:
    int Z0;
    int N0;
    int S0;
    int E0;
    type* hA;
    type* dA;
    type* A;
    HDV ~kFArray1D() { if(hA) free(hA); }
    HDV  kFArray1D(): Z0(0),S0(0),N0(0),E0(0)
    {  hA = NULL; dA = NULL; A = NULL; }
    HDV  kFArray1D(int n0): Z0(0),S0(0), N0(n0),E0(n0-1)
    {  hA = (type *)malloc(N0*sizeof(type)); dA = NULL; A = hA; }
    HDV  kFArray1D(int L0, int U0): Z0(0),S0(L0),N0(U0-L0+1),E0(U0)
    {  hA = (type *)malloc(N0*sizeof(type)); dA = NULL; A = hA; }

    HST kFArray1D<type>* clone()
    {   // copy structure without new allocation
        kFArray1D<type>* hM = new kFArray1D();
        hM->Z0 = Z0;
        hM->S0 = S0;    hM->E0 = E0;    hM->N0 = N0;
        hM->hA = hA;    hM->A  =  A;
        return hM;
    }
#if NVCC
    HST kFArray1D<type>* makeDeviceCopy()
    {
        int Nelmts = N0;
        type* saved_hA = hA;
        if(dA) cudaFree(dA); cudaMalloc((void **) &dA,Nelmts*sizeof(type));
        cudaMemcpy(dA, hA, Nelmts*sizeof(type), cudaMemcpyHostToDevice);
        hA = NULL; A = dA;
        kFArray1D<type>* dM;
        cudaMalloc((void **) &dM,sizeof(*this));
        cudaMemcpy(dM, this, sizeof(*this), cudaMemcpyHostToDevice);
        hA = saved_hA; A = hA;
        return dM;
    }
    HST kFArray1D<type>* makeDeviceOnly()
    {
        int Nelmts = N0;
        type* saved_hA = hA;
        if(dA) cudaFree(dA); cudaMalloc((void **) &dA,Nelmts*sizeof(type));
        // no copying of Array, only pointer info
        hA = NULL; A = dA;
        kFArray1D<type>* dM;
        cudaMalloc((void **) &dM,sizeof(*this));
        cudaMemcpy(dM, this, sizeof(*this), cudaMemcpyHostToDevice);
        hA = saved_hA; A = hA;
        return dM;
    }
    HST kFArray1D<type>* makeDeviceCopy(bool cpyArray)
    {
        if(cpyArray) return makeDeviceCopy();
        else         return makeDeviceOnly();
    }
    HST void freeDeviceCopy()  {  A=hA; if(dA) cudaFree(dA); dA = NULL; }
    HST void freeDeviceCopy(kFArray1D<type>* dM)  {  A=hA; if(dA) cudaFree(dA); dA = NULL; cudaFree(dM); }
    HST void copyDataHostToDevice() { cudaMemcpy(dA, hA, N0*sizeof(type), cudaMemcpyHostToDevice); }
    HST void copyDataDeviceToHost() { cudaMemcpy(hA, dA, N0*sizeof(type), cudaMemcpyDeviceToHost); }
#endif
    HDV inline type get(int i0)   {  return A[i0-S0]; }
    HDV inline void set(int i0, type v)  {  A[i0-S0] =v; }
    HDV inline void add(int i0, type v)  {  A[i0-S0]+=v; }
    HDV inline void mul(int i0, type v)  {  A[i0-S0]*=v; }
    HDV inline type&  operator()(int i0) {  return A[(i0-S0)]; }
    HDV inline void set(type  v)  {  for(int i=0;i<N0;i++) A[i] = v; }
    HDV inline void set(type* B)  {  for(int i=0;i<N0;i++) A[i] = B[i]; }
    HDV inline int size()   {  return N0; }
    HDV inline int length() {  return N0; }
    HDV inline int length(int d) {  return N0; }
    HDV inline type& at(int p) {  return A[p]; }
    HDV inline index1D indices(int p)
    {
        int Na = N0;
        int ia = p;
        index1D ndx;
        ia %= Na; Na /= N0; ndx.i0 = (int) (ia/Na) + S0;
        return ndx;
    }

    HDV void make(int n0)
    {   Z0= 0; S0= 0; N0= n0; E0= n0-1;
        if(hA) free(hA); hA = (type *)malloc(N0*sizeof(type));
        dA = NULL; A = hA;
    }
    HDV void make(int L0, int U0)
    {   Z0= 0; S0=L0; N0=U0-L0+1; E0=U0;
        if(hA) free(hA); hA = (type *)malloc(N0*sizeof(type));
        dA = NULL; A = hA;
    }

    HST void write(FILE* f)
    {
        fwrite(&Z0,sizeof(int),1,f);
        fwrite(&S0,sizeof(int),1,f);
        fwrite(&N0,sizeof(int),1,f);
        fwrite(hA,sizeof(type),N0,f);
    }
    HST void read (FILE* f)
    {   if(hA) free(hA);
        fread (&Z0,sizeof(int),1,f);
        fread (&S0,sizeof(int),1,f);
        fread (&N0,sizeof(int),1,f);
        hA = (type *)malloc(N0*sizeof(type));
        fread(hA,sizeof(type),N0,f);
        A = hA;
        E0=S0+N0-1;
    }
};


template<class type> class kFArray2D
{ // col-major storage: first index varies first
public:
    int Z0;
    int N0,N1;  // array size in each direction
    int S0,S1;  // index shift in each direction
    int E0,E1;  // maximum-index in each direction
    type* hA;
    type* dA;
    type* A;
    HDV  ~kFArray2D() { if(hA) free(hA); }
    HDV   kFArray2D(): Z0(0),S0(0),S1(0),
                             N0(0),N1(0),
                             E0(0),E1(0)
    { hA = NULL; dA = NULL; A = NULL; }
    HDV   kFArray2D(int n0, int n1):Z0(0),
                                    S0(0),   S1(0),
                                    N0(n0),  N1(n1),
                                    E0(n0-1),E1(n1-1)
    {  hA = (type *)malloc(N0*N1*sizeof(type)); dA = NULL; A = hA; }
    HDV   kFArray2D(int L0, int U0,
                    int L1, int U1):Z0(0),
                                    S0(L0),     S1(L1),
                                    N0(U0-L0+1),N1(U1-L1+1),
                                    E0(U0),     E1(U1)
    {  hA = (type *)malloc(N0*N1*sizeof(type)); dA = NULL; A = hA; }

    HST kFArray2D<type>* clone()
    {   // copy structure without new allocation
        kFArray2D<type>* hM = new kFArray2D();
        hM->Z0 = Z0;
        hM->S0 = S0;    hM->E0 = E0;    hM->N0 = N0;
        hM->S1 = S1;    hM->E1 = E1;    hM->N1 = N1;
        hM->hA = hA;    hM->A  =  A;
        return hM;
    }
#if NVCC
    HST kFArray2D<type>* makeDeviceCopy()
    {
        int Nelmts = N0*N1;
        type* saved_hA = hA;
        if(dA) cudaFree(dA); cudaMalloc((void **) &dA,Nelmts*sizeof(type));
        cudaMemcpy(dA, hA, Nelmts*sizeof(type), cudaMemcpyHostToDevice);
        hA = NULL; A = dA;
        kFArray2D<type>* dM;
        cudaMalloc((void **) &dM,sizeof(*this));
        cudaMemcpy(dM, this, sizeof(*this), cudaMemcpyHostToDevice);
        hA = saved_hA; A = hA;
        return dM;
    }
    HST kFArray2D<type>* makeDeviceOnly()
    {
        int Nelmts = N0*N1;
        type* saved_hA = hA;
        if(dA) cudaFree(dA); cudaMalloc((void **) &dA,Nelmts*sizeof(type));
        // no copying of Array, only pointer info
        hA = NULL; A = dA;
        kFArray2D<type>* dM;
        cudaMalloc((void **) &dM,sizeof(*this));
        cudaMemcpy(dM, this, sizeof(*this), cudaMemcpyHostToDevice);
        hA = saved_hA; A = hA;
        return dM;
    }
    HST kFArray2D<type>* makeDeviceCopy(bool cpyArray)
    {
        if(cpyArray) return makeDeviceCopy();
        else         return makeDeviceOnly();
    }
    HST void freeDeviceCopy()  {  A=hA; if(dA) cudaFree(dA); dA = NULL; }
    HST void freeDeviceCopy(kFArray2D<type>* dM)  {  A=hA; if(dA) cudaFree(dA); dA = NULL; cudaFree(dM); }
    HST void copyDataHostToDevice() { cudaMemcpy(dA, hA, N0*N1*sizeof(type), cudaMemcpyHostToDevice); }
    HST void copyDataDeviceToHost() { cudaMemcpy(hA, dA, N0*N1*sizeof(type), cudaMemcpyDeviceToHost); }
#endif
    HDV inline type get(int i0, int i1)   {  return A[(i1-S1)*N0+i0-S0]; }
    HDV inline void set(int i0, int i1, type v)  {  A[(i1-S1)*N0+i0-S0] =v; }
    HDV inline void add(int i0, int i1, type v)  {  A[(i1-S1)*N0+i0-S0]+=v; }
    HDV inline void mul(int i0, int i1, type v)  {  A[(i1-S1)*N0+i0-S0]*=v; }
    HDV inline type&  operator()(int i0, int i1) {  return A[(i1-S1)*N0+i0-S0]; }
    HDV inline void set(type  v)  {  for(int i=0;i<N1*N0;i++) A[i] = v; }
    HDV inline void set(type *B)  {  for(int i=0;i<N1*N0;i++) A[i] = B[i]; }
    HDV inline int size() {  return N0*N1; }
    HDV inline int length(int d) {  if(d==0) return N0; if(d==1) return N1; }
    HDV inline type& at(int p) {  return A[p]; }
    HDV inline index2D indices(int p)
    {
        int Na = N1*N0;
        int ia = p;
        index2D ndx;
        ia %= Na; Na /= N1; ndx.i1 = (int) (ia/Na) + S1;
        ia %= Na; Na /= N0; ndx.i0 = (int) (ia/Na) + S0;
        return ndx;
    }

    HDV kFArray1D<type>& subMlast(int i1);  // make 1D array from last  index i1
    HDV kFArray1D<type>& subMfrst(int i0);  // make 1D array from first index i0
    HDV kFArray1D<type>* psubMlast(int i1);  // make 1D array from last  index i1
    HDV kFArray1D<type>* psubMfrst(int i0);  // make 1D array from first index i0

    HDV void make(int n0, int n1)
    {   Z0= 0; S0= 0; N0= n0; E0= n0-1;
        S1= 0; N1= n1; E1= n1-1;
        if(hA) free(hA); hA = (type *)malloc(N1*N0*sizeof(type));
        dA = NULL; A = hA;
    }
    HDV void make(int L0, int U0,
                                  int L1, int U1)
    {   Z0= 0; S0=L0; N0=U0-L0+1; E0=U0;
        S1=L1; N1=U1-L1+1; E1=U1;
        if(hA) free(hA); hA = (type *)malloc(N1*N0*sizeof(type));
        dA = NULL; A = hA;
    }

    HST void write(FILE* f)
    {
        fwrite(&Z0,sizeof(int),1,f);
        fwrite(&S0,sizeof(int),1,f);
        fwrite(&N0,sizeof(int),1,f);
        fwrite(&S1,sizeof(int),1,f);
        fwrite(&N1,sizeof(int),1,f);
        fwrite(hA,sizeof(type),N0*N1,f);
    }
    HST void read (FILE* f)
    {   if(hA) free(hA);
        fread (&Z0,sizeof(int),1,f);
        fread (&S0,sizeof(int),1,f);
        fread (&N0,sizeof(int),1,f);
        fread (&S1,sizeof(int),1,f);
        fread (&N1,sizeof(int),1,f);
        hA = (type *)malloc(N0*N1*sizeof(type));
        fread(hA,sizeof(type),N0*N1,f);
        A = hA;
        E0=S0+N0-1;
        E1=S1+N1-1;
    }
};


template<class type> class kFArray3D
{ // col-major storage: first index varies first
public:
    int Z0;
    int N0,N1,N2;
    int S0,S1,S2;
    int E0,E1,E2;
    type* hA;
    type* dA;
    type* A;
    HDV ~kFArray3D() { if(hA) free(hA); }
    HDV  kFArray3D():Z0(0),
                     S0(0),S1(0),S2(0),
                     N0(0),N1(0),N2(0),
                     E0(0),E1(0),E2(0)
    {  hA = NULL; dA = NULL; A = NULL;}
    HDV  kFArray3D(int n0, int n1, int n2): Z0(0),
                                   S0(0),   S1(0),   S2(0),
                                   N0(n0),  N1(n1),  N2(n2),
                                   E0(n0-1),E1(n1-1),E2(n2-1)
    {  hA = (type *)malloc(N0*N1*N2*sizeof(type)); dA = NULL; A = hA; }
    HDV  kFArray3D(int L0, int U0,
                   int L1, int U1,
                   int L2, int U2):Z0(0),
                                   S0(L0),     S1(L1),     S2(L2),
                                   N0(U0-L0+1),N1(U1-L1+1),N2(U2-L2+1),
                                   E0(U0),     E1(U1),     E2(U2)
    {  hA = (type *)malloc(N0*N1*N2*sizeof(type)); dA = NULL; A = hA; }

    HST kFArray3D<type>* clone()
    {   // copy structure without new allocation
        kFArray3D<type>* hM = new kFArray3D();
        hM->Z0 = Z0;
        hM->S0 = S0;    hM->E0 = E0;    hM->N0 = N0;
        hM->S1 = S1;    hM->E1 = E1;    hM->N1 = N1;
        hM->S2 = S2;    hM->E2 = E2;    hM->N2 = N2;
        hM->hA = hA;    hM->A  =  A;
        return hM;
    }

#if NVCC
    HST kFArray3D<type>* makeDeviceCopy()
    {
        int Nelmts = N0*N1*N2;
        type* saved_hA = hA;
        if(dA) cudaFree(dA); cudaMalloc((void **) &dA,Nelmts*sizeof(type));
        cudaMemcpy(dA, hA, Nelmts*sizeof(type), cudaMemcpyHostToDevice);
        hA = NULL; A = dA;
        kFArray3D<type>* dM;
        cudaMalloc((void **) &dM,sizeof(*this));
        cudaMemcpy(dM, this, sizeof(*this), cudaMemcpyHostToDevice);
        hA = saved_hA; A = hA;
        return dM;
    }
    HST kFArray3D<type>* makeDeviceOnly()
    {
        int Nelmts = N0*N1*N2;
        type* saved_hA = hA;
        if(dA) cudaFree(dA); cudaMalloc((void **) &dA,Nelmts*sizeof(type));
        // no copying of Array, only pointer info
        hA = NULL; A = dA;
        kFArray3D<type>* dM;
        cudaMalloc((void **) &dM,sizeof(*this));
        cudaMemcpy(dM, this, sizeof(*this), cudaMemcpyHostToDevice);
        hA = saved_hA; A = hA;
        return dM;
    }
    HST kFArray3D<type>* makeDeviceCopy(bool cpyArray)
    {
        if(cpyArray) return makeDeviceCopy();
        else         return makeDeviceOnly();
    }
    HST void freeDeviceCopy()  {  A=hA; if(dA) cudaFree(dA); dA = NULL; }
    HST void freeDeviceCopy(kFArray3D<type>* dM)  {  A=hA; if(dA) cudaFree(dA); dA = NULL; cudaFree(dM); }
    HST void copyDataHostToDevice() { cudaMemcpy(dA, hA, N0*N1*N2*sizeof(type), cudaMemcpyHostToDevice); }
    HST void copyDataDeviceToHost() { cudaMemcpy(hA, dA, N0*N1*N2*sizeof(type), cudaMemcpyDeviceToHost); }
#endif
    HDV inline type get(int i0, int i1, int i2)   {  return A[((i2-S2)*N1+i1-S1)*N0+i0-S0]; }
    HDV inline void set(int i0, int i1, int i2, type v)  {  A[((i2-S2)*N1+i1-S1)*N0+i0-S0] =v; }
    HDV inline void add(int i0, int i1, int i2, type v)  {  A[((i2-S2)*N1+i1-S1)*N0+i0-S0]+=v; }
    HDV inline void mul(int i0, int i1, int i2, type v)  {  A[((i2-S2)*N1+i1-S1)*N0+i0-S0]*=v; }
    HDV inline type&  operator()(int i0, int i1, int i2) { return A[((i2-S2)*N1+i1-S1)*N0+i0-S0]; }
    HDV inline void set(type  v)  {  for(int i=0;i<N2*N1*N0;i++) A[i] = v; }
    HDV inline void set(type *B)  {  for(int i=0;i<N2*N1*N0;i++) A[i] = B[i]; }
    HDV inline int size() {  return N0*N1*N2; }
    HDV inline int length(int d) {  if(d==0) return N0; if(d==1) return N1; if(d==2) return N2; }
    HDV inline type& at(int p) {  return A[p]; }
    HDV inline index3D indices(int p)
    {
        int Na = N2*N1*N0;
        int ia = p;
        index3D ndx;
        ia %= Na; Na /= N2; ndx.i2 = (int) (ia/Na) + S2;
        ia %= Na; Na /= N1; ndx.i1 = (int) (ia/Na) + S1;
        ia %= Na; Na /= N0; ndx.i0 = (int) (ia/Na) + S0;
        return ndx;
    }

    HDV kFArray2D<type>& subMlast(int i2);  // make 2D array from last  index i2
    HDV kFArray2D<type>& subMfrst(int i0);  // make 2D array from first index i0

    HDV void make(int n0, int n1, int n2)
    {   Z0= 0;
        S0= 0; N0= n0; E0= n0-1;
        S1= 0; N1= n1; E1= n1-1;
        S2= 0; N2= n2; E2= n2-1;
        if(hA) free(hA); hA = (type *)malloc(N2*N1*N0*sizeof(type));
        dA = NULL; A = hA;
    }
    HDV void make(int L0, int U0,
                                  int L1, int U1,
                                  int L2, int U2)
    {   Z0= 0;
        S0=L0; N0=U0-L0+1; E0=U0;
        S1=L1; N1=U1-L1+1; E1=U1;
        S2=L2; N2=U2-L2+1; E2=U2;
        if(hA) free(hA); hA = (type *)malloc(N2*N1*N0*sizeof(type));
        dA = NULL; A = hA;
    }

    HST void write(FILE* f)
    {
        fwrite(&Z0,sizeof(int),1,f);
        fwrite(&S0,sizeof(int),1,f);
        fwrite(&N0,sizeof(int),1,f);
        fwrite(&S1,sizeof(int),1,f);
        fwrite(&N1,sizeof(int),1,f);
        fwrite(&S2,sizeof(int),1,f);
        fwrite(&N2,sizeof(int),1,f);
        fwrite(hA,sizeof(type),N0*N1*N2,f);
    }
    HST void read (FILE* f)
    {   if(hA) free(hA);
        fread (&Z0,sizeof(int),1,f);
        fread (&S0,sizeof(int),1,f);
        fread (&N0,sizeof(int),1,f);
        fread (&S1,sizeof(int),1,f);
        fread (&N1,sizeof(int),1,f);
        fread (&S2,sizeof(int),1,f);
        fread (&N2,sizeof(int),1,f);
        hA = (type *)malloc(N0*N1*N2*sizeof(type));
        fread(hA,sizeof(type),N0*N1*N2,f);
        A = hA;
        E0=S0+N0-1;
        E1=S1+N1-1;
        E2=S2+N2-1;
    }
};

template<class type> class kFArray4D
{ // col-major storage: first index varies first
public:
    int Z0;
    int N0,N1,N2,N3;
    int S0,S1,S2,S3;
    int E0,E1,E2,E3;
    type* hA;
    type* dA;
    type* A;
    HDV  ~kFArray4D() { if(hA) free(hA); }
    HDV   kFArray4D():Z0(0),
                      S0(0),S1(0),S2(0),S3(0),
                      N0(0),N1(0),N2(0),N3(0),
                      E0(0),E1(0),E2(0),E3(0)
    { hA = NULL; dA = NULL; A = NULL; }
    HDV   kFArray4D(int n0, int n1, int n2, int n3): Z0(0),
                                    S0( 0)  ,S1( 0)  ,S2( 0)  ,S3( 0)  ,
                                    N0(n0  ),N1(n1  ),N2(n2  ),N3(n3  ),
                                    E0(n0-1),E1(n1-1),E2(n2-1),E3(n3-1)
    {  hA = (type *)malloc(N0*N1*N2*N3*sizeof(type)); dA = NULL; A = hA; }
    HDV   kFArray4D(int L0, int U0,
                    int L1, int U1,
                    int L2, int U2,
                    int L3, int U3):Z0(0),
                                    S0(L0)     ,S1(L1)     ,S2(L2)     ,S3(L3),
                                    N0(U0-L0+1),N1(U1-L1+1),N2(U2-L2+1),N3(U3-L3+1),
                                    E0(U0)     ,E1(U1)     ,E2(U2)     ,E3(U3)
    {  hA = (type *)malloc(N0*N1*N2*N3*sizeof(type)); dA = NULL; A = hA; }

    HST kFArray4D<type>* clone()
    {   // copy structure without new allocation
        kFArray4D<type>* hM = new kFArray4D();
        hM->Z0 = Z0;
        hM->S0 = S0;    hM->E0 = E0;    hM->N0 = N0;
        hM->S1 = S1;    hM->E1 = E1;    hM->N1 = N1;
        hM->S2 = S2;    hM->E2 = E2;    hM->N2 = N2;
        hM->S3 = S3;    hM->E3 = E3;    hM->N3 = N3;
        hM->hA = hA;    hM->A  =  A;
        return hM;
    }
#if NVCC
    HST kFArray4D<type>* makeDeviceCopy()
    {
        int Nelmts = N0*N1*N2*N3;
        type* saved_hA = hA;
        if(dA) cudaFree(dA); cudaMalloc((void **) &dA,Nelmts*sizeof(type));
        cudaMemcpy(dA, hA, Nelmts*sizeof(type), cudaMemcpyHostToDevice);
        hA = NULL; A = dA;
        kFArray4D<type>* dM;
        cudaMalloc((void **) &dM,sizeof(*this));
        cudaMemcpy(dM, this, sizeof(*this), cudaMemcpyHostToDevice);
        hA = saved_hA; A = hA;
        return dM;
    }
    HST kFArray4D<type>* makeDeviceOnly()
    {
        int Nelmts = N0*N1*N2*N3;
        type* saved_hA = hA;
        if(dA) cudaFree(dA); cudaMalloc((void **) &dA,Nelmts*sizeof(type));
        // no copying of Array, only pointer info
        hA = NULL; A = dA;
        kFArray4D<type>* dM;
        cudaMalloc((void **) &dM,sizeof(*this));
        cudaMemcpy(dM, this, sizeof(*this), cudaMemcpyHostToDevice);
        hA = saved_hA; A = hA;
        return dM;
    }
    HST kFArray4D<type>* makeDeviceCopy(bool cpyArray)
    {
        if(cpyArray) return makeDeviceCopy();
        else         return makeDeviceOnly();
    }
    HST void freeDeviceCopy()  {  A=hA; if(dA) cudaFree(dA); dA = NULL; }
    HST void freeDeviceCopy(kFArray4D<type>* dM)  {  A=hA; if(dA) cudaFree(dA); dA = NULL; cudaFree(dM); }
    HST void copyDataHostToDevice() { cudaMemcpy(dA, hA, N0*N1*N2*N3*sizeof(type), cudaMemcpyHostToDevice); }
    HST void copyDataDeviceToHost() { cudaMemcpy(hA, dA, N0*N1*N2*N3*sizeof(type), cudaMemcpyDeviceToHost); }
#endif
    HDV inline type get(int i0, int i1, int i2, int i3)   {  return A[(((i3-S3)*N2+i2-S2)*N1+i1-S1)*N0+i0-S0]; }
    HDV inline void set(int i0, int i1, int i2, int i3, type v)  {  A[(((i3-S3)*N2+i2-S2)*N1+i1-S1)*N0+i0-S0] =v; }
    HDV inline void add(int i0, int i1, int i2, int i3, type v)  {  A[(((i3-S3)*N2+i2-S2)*N1+i1-S1)*N0+i0-S0]+=v; }
    HDV inline void mul(int i0, int i1, int i2, int i3, type v)  {  A[(((i3-S3)*N2+i2-S2)*N1+i1-S1)*N0+i0-S0]*=v; }
    HDV inline type&  operator()(int i0, int i1, int i2, int i3) {  return A[(((i3-S3)*N2+i2-S2)*N1+i1-S1)*N0+i0-S0]; }
    HDV inline void set(type  v)  {  for(int i=0;i<N3*N2*N1*N0;i++) A[i] = v; }
    HDV inline void set(type *B)  {  for(int i=0;i<N3*N2*N1*N0;i++) A[i] = B[i]; }
    HDV inline int size() {  return N0*N1*N2*N3; }
    HDV inline int length(int d) {  if(d==0) return N0; if(d==1) return N1; if(d==2) return N2; if(d==3) return N3; }
    HDV inline type& at(int p) {  return A[p]; }
    HDV inline index4D indices(int p)
    {
        int Na = N3*N2*N1*N0;
        int ia = p;
        index4D ndx;
        ia %= Na; Na /= N3; ndx.i3 = (int) (ia/Na) + S3;
        ia %= Na; Na /= N2; ndx.i2 = (int) (ia/Na) + S2;
        ia %= Na; Na /= N1; ndx.i1 = (int) (ia/Na) + S1;
        ia %= Na; Na /= N0; ndx.i0 = (int) (ia/Na) + S0;
        return ndx;
    }

    HDV kFArray3D<type>& subMlast(int i3);  // make 3D array from last  index i3
    HDV kFArray3D<type>& subMfrst(int i0);  // make 3D array from first index i0

    HDV void make(int n0, int n1, int n2, int n3)
    {   Z0= 0;
        S0= 0; N0= n0; E0= n0-1;
        S1= 0; N1= n1; E1= n1-1;
        S2= 0; N2= n2; E2= n2-1;
        S3= 0; N3= n3; E3= n3-1;
        if(hA) free(hA); hA = (type *)malloc(N3*N2*N1*N0*sizeof(type));
        dA = NULL; A = hA;
    }
    HDV void make(int L0, int U0,
                                  int L1, int U1,
                                  int L2, int U2,
                                  int L3, int U3)
    {   Z0= 0;
        S0=L0; N0=U0-L0+1; E0=U0;
        S1=L1; N1=U1-L1+1; E1=U1;
        S2=L2; N2=U2-L2+1; E2=U2;
        S3=L3; N3=U3-L3+1; E3=U3;
        if(hA) free(hA); hA = (type *)malloc(N3*N2*N1*N0*sizeof(type));
        dA = NULL; A = hA;
    }

    HST void write(FILE* f)
    {
        fwrite(&Z0,sizeof(int),1,f);
        fwrite(&S0,sizeof(int),1,f);
        fwrite(&N0,sizeof(int),1,f);
        fwrite(&S1,sizeof(int),1,f);
        fwrite(&N1,sizeof(int),1,f);
        fwrite(&S2,sizeof(int),1,f);
        fwrite(&N2,sizeof(int),1,f);
        fwrite(&S3,sizeof(int),1,f);
        fwrite(&N3,sizeof(int),1,f);
        fwrite(hA,sizeof(type),N0*N1*N2*N3,f);
    }
    HST void read (FILE* f)
    {   if(hA) free(hA);
        fread (&Z0,sizeof(int),1,f);
        fread (&S0,sizeof(int),1,f);
        fread (&N0,sizeof(int),1,f);
        fread (&S1,sizeof(int),1,f);
        fread (&N1,sizeof(int),1,f);
        fread (&S2,sizeof(int),1,f);
        fread (&N2,sizeof(int),1,f);
        fread (&S3,sizeof(int),1,f);
        fread (&N3,sizeof(int),1,f);
        hA = (type *)malloc(N0*N1*N2*N3*sizeof(type));
        fread(hA,sizeof(type),N0*N1*N2*N3,f);
        A = hA;
        E0=S0+N0-1;
        E1=S1+N1-1;
        E2=S2+N2-1;
        E3=S3+N3-1;
    }
};

template<class type> class kFArray5D
{ // col-major storage: first index varies first
public:
    int Z0;
    int N0,N1,N2,N3,N4;
    int S0,S1,S2,S3,S4;
    int E0,E1,E2,E3,E4;
    type* hA;
    type* dA;
    type* A;
    HDV ~kFArray5D() { if(hA) free(hA); }
    HDV  kFArray5D():Z0(0),
                     S0(0),S1(0),S2(0),S3(0),S4(0),
                     N0(0),N1(0),N2(0),N3(0),N4(0),
                     E0(0),E1(0),E2(0),E3(0),E4(0)
    {  hA = NULL; dA = NULL;A = NULL; }
    HDV  kFArray5D(int n0, int n1, int n2, int n3, int n4): Z0(0),
                                   S0(0)   ,S1(0)   ,S2(0)   ,S3(0)   ,S4(0),
                                   N0(n0  ),N1(n1  ),N2(n2  ),N3(n3  ),N4(n4  ),
                                   E0(n0-1),E1(n1-1),E2(n2-1),E3(n3-1),E4(n4-1)
    {  hA = (type *)malloc(N0*N1*N2*N3*N4*sizeof(type)); dA = NULL; A = hA; }
    HDV  kFArray5D(int L0, int U0,
                   int L1, int U1,
                   int L2, int U2,
                   int L3, int U3,
                   int L4, int U4):Z0(0),
                                   S0(L0)     ,S1(L1)     ,S2(L2)     ,S3(L3)     ,S4(L4),
                                   N0(U0-L0+1),N1(U1-L1+1),N2(U2-L2+1),N3(U3-L3+1),N4(U4-L4+1),
                                   E0(U0)     ,E1(U1)     ,E2(U2)     ,E3(U3)     ,E4(U4)
    {  hA = (type *)malloc(N0*N1*N2*N3*N4*sizeof(type)); dA = NULL; A = hA; }

    HST kFArray5D<type>* clone()
    {   // copy structure without new allocation
        kFArray5D<type>* hM = new kFArray5D();
        hM->Z0 = Z0;
        hM->S0 = S0;    hM->E0 = E0;    hM->N0 = N0;
        hM->S1 = S1;    hM->E1 = E1;    hM->N1 = N1;
        hM->S2 = S2;    hM->E2 = E2;    hM->N2 = N2;
        hM->S3 = S3;    hM->E3 = E3;    hM->N3 = N3;
        hM->S4 = S4;    hM->E4 = E4;    hM->N4 = N4;
        hM->hA = hA;    hM->A  =  A;
        return hM;
    }
#if NVCC
    HST kFArray5D<type>* makeDeviceCopy()
    {
        int Nelmts = N0*N1*N2*N3*N4;
        type* saved_hA = hA;
        if(dA) cudaFree(dA); cudaMalloc((void **) &dA,Nelmts*sizeof(type));
        cudaMemcpy(dA, hA, Nelmts*sizeof(type), cudaMemcpyHostToDevice);
        hA = NULL; A = dA;
        kFArray5D<type>* dM;
        cudaMalloc((void **) &dM,sizeof(*this));
        cudaMemcpy(dM, this, sizeof(*this), cudaMemcpyHostToDevice);
        hA = saved_hA; A = hA;
        return dM;
    }
    HST kFArray5D<type>* makeDeviceOnly()
    {
        int Nelmts = N0*N1*N2*N3*N4;
        type* saved_hA = hA;
        if(dA) cudaFree(dA); cudaMalloc((void **) &dA,Nelmts*sizeof(type));
        // no copying of Array, only pointer info
        hA = NULL; A = dA;
        kFArray5D<type>* dM;
        cudaMalloc((void **) &dM,sizeof(*this));
        cudaMemcpy(dM, this, sizeof(*this), cudaMemcpyHostToDevice);
        hA = saved_hA; A = hA;
        return dM;
    }
    HST kFArray5D<type>* makeDeviceCopy(bool cpyArray)
    {
        if(cpyArray) return makeDeviceCopy();
        else         return makeDeviceOnly();
    }

    HST void freeDeviceCopy()  {  A=hA; if(dA) cudaFree(dA); dA = NULL; }
    HST void freeDeviceCopy(kFArray5D<type>* dM)  {  A=hA; if(dA) cudaFree(dA); dA = NULL; cudaFree(dM); }
    HST void copyDataHostToDevice() { cudaMemcpy(dA, hA, N0*N1*N2*N3*N4*sizeof(type), cudaMemcpyHostToDevice); }
    HST void copyDataDeviceToHost() { cudaMemcpy(hA, dA, N0*N1*N2*N3*N4*sizeof(type), cudaMemcpyDeviceToHost); }
#endif
    HDV inline type get(int i0, int i1, int i2, int i3, int i4)   {  return A[((((i4-S4)*N3+i3-S3)*N2+i2-S2)*N1+i1-S1)*N0+i0-S0]; }
    HDV inline void set(int i0, int i1, int i2, int i3, int i4, type v)  {  A[((((i4-S4)*N3+i3-S3)*N2+i2-S2)*N1+i1-S1)*N0+i0-S0] =v; }
    HDV inline void add(int i0, int i1, int i2, int i3, int i4, type v)  {  A[((((i4-S4)*N3+i3-S3)*N2+i2-S2)*N1+i1-S1)*N0+i0-S0]+=v; }
    HDV inline void mul(int i0, int i1, int i2, int i3, int i4, type v)  {  A[((((i4-S4)*N3+i3-S3)*N2+i2-S2)*N1+i1-S1)*N0+i0-S0]*=v; }
    HDV inline type&  operator()(int i0, int i1, int i2, int i3, int i4) {  return A[((((i4-S4)*N3+i3-S3)*N2+i2-S2)*N1+i1-S1)*N0+i0-S0]; }
    HDV inline void set(type  v)  {  for(int i=0;i<N4*N3*N2*N1*N0;i++) A[i] = v; }
    HDV inline void set(type *B)  {  for(int i=0;i<N4*N3*N2*N1*N0;i++) A[i] = B[i]; }
    HDV inline int size() {  return N0*N1*N2*N3*N4; }
    HDV inline int length(int d) {  if(d==0) return N0; if(d==1) return N1; if(d==2) return N2; if(d==3) return N3;
                                    if(d==4) return N4; }
    HDV inline type& at(int p) {  return A[p]; }
    HDV inline index5D indices(int p)
    {
        int Na = N4*N3*N2*N1*N0;
        int ia = p;
        index5D ndx;
        ia %= Na; Na /= N4; ndx.i4 = (int) (ia/Na) + S4;
        ia %= Na; Na /= N3; ndx.i3 = (int) (ia/Na) + S3;
        ia %= Na; Na /= N2; ndx.i2 = (int) (ia/Na) + S2;
        ia %= Na; Na /= N1; ndx.i1 = (int) (ia/Na) + S1;
        ia %= Na; Na /= N0; ndx.i0 = (int) (ia/Na) + S0;
        return ndx;
    }

    HDV kFArray4D<type>& subMlast(int i4);  // make 4D array from last  index i4
    HDV kFArray4D<type>& subMfrst(int i0);  // make 4D array from first index i0

    HDV void make(int n0, int n1, int n2, int n3, int n4)
    {   Z0= 0;
        S0= 0; N0= n0; E0= n0-1;
        S1= 0; N1= n1; E1= n1-1;
        S2= 0; N2= n2; E2= n2-1;
        S3= 0; N3= n3; E3= n3-1;
        S4= 0; N4= n4; E4= n4-1;
        if(hA) free(hA); hA = (type *)malloc(N4*N3*N2*N1*N0*sizeof(type));
        dA = NULL; A = hA;
    }
    HDV void make(int L0, int U0,
                                  int L1, int U1,
                                  int L2, int U2,
                                  int L3, int U3,
                                  int L4, int U4)
    {   Z0= 0;
        S0=L0; N0=U0-L0+1; E0=U0;
        S1=L1; N1=U1-L1+1; E1=U1;
        S2=L2; N2=U2-L2+1; E2=U2;
        S3=L3; N3=U3-L3+1; E3=U3;
        S4=L4; N4=U4-L4+1; E4=U4;
        if(hA) free(hA); hA = (type *)malloc(N4*N3*N2*N1*N0*sizeof(type));
        dA = NULL; A = hA;
    }

    HST void write(FILE* f)
    {
        fwrite(&Z0,sizeof(int),1,f);
        fwrite(&S0,sizeof(int),1,f);
        fwrite(&N0,sizeof(int),1,f);
        fwrite(&S1,sizeof(int),1,f);
        fwrite(&N1,sizeof(int),1,f);
        fwrite(&S2,sizeof(int),1,f);
        fwrite(&N2,sizeof(int),1,f);
        fwrite(&S3,sizeof(int),1,f);
        fwrite(&N3,sizeof(int),1,f);
        fwrite(&S4,sizeof(int),1,f);
        fwrite(&N4,sizeof(int),1,f);
        fwrite(hA,sizeof(type),N0*N1*N2*N3*N4,f);
    }
    HST void read (FILE* f)
    {   if(hA) free(hA);
        fread (&Z0,sizeof(int),1,f);
        fread (&S0,sizeof(int),1,f);
        fread (&N0,sizeof(int),1,f);
        fread (&S1,sizeof(int),1,f);
        fread (&N1,sizeof(int),1,f);
        fread (&S2,sizeof(int),1,f);
        fread (&N2,sizeof(int),1,f);
        fread (&S3,sizeof(int),1,f);
        fread (&N3,sizeof(int),1,f);
        fread (&S4,sizeof(int),1,f);
        fread (&N4,sizeof(int),1,f);
        hA = (type *)malloc(N0*N1*N2*N3*N4*sizeof(type));
        fread(hA,sizeof(type),N0*N1*N2*N3*N4,f);
        A = hA;
        E0=S0+N0-1;
        E1=S1+N1-1;
        E2=S2+N2-1;
        E3=S3+N3-1;
        E4=S4+N4-1;
    }
};

template<class type> class kFArray6D
{ // col-major storage: first index varies first
public:
    int Z0;
    int N0,N1,N2,N3,N4,N5;
    int S0,S1,S2,S3,S4,S5;
    int E0,E1,E2,E3,E4,E5;
    type* hA;
    type* dA;
    type* A;
    HDV ~kFArray6D() { if(hA) free(hA); }
    HDV  kFArray6D():Z0(0),
                     S0(0),S1(0),S2(0),S3(0),S4(0),S5(0),
                     N0(0),N1(0),N2(0),N3(0),N4(0),N5(0),
                     E0(0),E1(0),E2(0),E3(0),E4(0),E5(0)
    {  hA = NULL;dA = NULL; A = NULL; }
    HDV  kFArray6D(int n0, int n1, int n2, int n3, int n4, int n5): Z0(0),
                                   S0(0)   ,S1(0)   ,S2(0)   ,S3(0)   ,S4(0)   ,S5(0),
                                   N0(n0  ),N1(n1  ),N2(n2  ),N3(n3  ),N4(n4  ),N5(n5  ),
                                   E0(n0-1),E1(n1-1),E2(n2-1),E3(n3-1),E4(n4-1),E5(n5-1)
    {  hA = (type *)malloc(N0*N1*N2*N3*N4*N5*sizeof(type)); dA = NULL; A = hA; }
    HDV  kFArray6D(int L0, int U0,
                   int L1, int U1,
                   int L2, int U2,
                   int L3, int U3,
                   int L4, int U4,
                   int L5, int U5):Z0(0),
                                   S0(L0)     ,S1(L1)     ,S2(L2)     ,S3(L3)     ,S4(L4)     ,S5(L5),
                                   N0(U0-L0+1),N1(U1-L1+1),N2(U2-L2+1),N3(U3-L3+1),N4(U4-L4+1),N5(U5-L5+1),
                                   E0(U0)     ,E1(U1)     ,E2(U2)     ,E3(U3)     ,E4(U4)     ,E5(U5)
    {  hA = (type *)malloc(N0*N1*N2*N3*N4*N5*sizeof(type)); dA = NULL; A = hA; }

    HST kFArray6D<type>* clone()
    {   // copy structure without new allocation
        kFArray6D<type>* hM = new kFArray6D();
        hM->Z0 = Z0;
        hM->S0 = S0;    hM->E0 = E0;    hM->N0 = N0;
        hM->S1 = S1;    hM->E1 = E1;    hM->N1 = N1;
        hM->S2 = S2;    hM->E2 = E2;    hM->N2 = N2;
        hM->S3 = S3;    hM->E3 = E3;    hM->N3 = N3;
        hM->S4 = S4;    hM->E4 = E4;    hM->N4 = N4;
        hM->S5 = S5;    hM->E5 = E5;    hM->N5 = N5;
        hM->hA = hA;    hM->A  =  A;
        return hM;
    }
#if NVCC
    HST kFArray6D<type>* makeDeviceCopy()
    {
        int Nelmts = N0*N1*N2*N3*N4*N5;
        type* saved_hA = hA;
        if(dA) cudaFree(dA); cudaMalloc((void **) &dA,Nelmts*sizeof(type));
        cudaMemcpy(dA, hA, Nelmts*sizeof(type), cudaMemcpyHostToDevice);
        hA = NULL; A = dA;
        kFArray6D<type>* dM;
        cudaMalloc((void **) &dM,sizeof(*this));
        cudaMemcpy(dM, this, sizeof(*this), cudaMemcpyHostToDevice);
        hA = saved_hA; A = hA;
        return dM;
    }
    HST kFArray6D<type>* makeDeviceOnly()
    {
        int Nelmts = N0*N1*N2*N3*N4*N5;
        type* saved_hA = hA;
        if(dA) cudaFree(dA); cudaMalloc((void **) &dA,Nelmts*sizeof(type));
        // no copying of Array, only pointer info
        hA = NULL; A = dA;
        kFArray6D<type>* dM;
        cudaMalloc((void **) &dM,sizeof(*this));
        cudaMemcpy(dM, this, sizeof(*this), cudaMemcpyHostToDevice);
        hA = saved_hA; A = hA;
        return dM;
    }

    HST kFArray6D<type>* makeDeviceCopy(bool cpyArray)
    {
        if(cpyArray) return makeDeviceCopy();
        else         return makeDeviceOnly();
    }
    HST void freeDeviceCopy()  {  A=hA; if(dA) cudaFree(dA); dA = NULL; }
    HST void freeDeviceCopy(kFArray6D<type>* dM)  {  A=hA; if(dA) cudaFree(dA); dA = NULL; cudaFree(dM); }
    HST void copyDataHostToDevice() { cudaMemcpy(dA, hA, N0*N1*N2*N3*N4*N5*sizeof(type), cudaMemcpyHostToDevice); }
    HST void copyDataDeviceToHost() { cudaMemcpy(hA, dA, N0*N1*N2*N3*N4*N5*sizeof(type), cudaMemcpyDeviceToHost); }
#endif
    HDV inline type get(int i0, int i1, int i2, int i3, int i4, int i5)   {  return A[(((((i5-S5)*N4+i4-S4)*N3+i3-S3)*N2+i2-S2)*N1+i1-S1)*N0+i0-S0]; }
    HDV inline void set(int i0, int i1, int i2, int i3, int i4, int i5, type v)  {  A[(((((i5-S5)*N4+i4-S4)*N3+i3-S3)*N2+i2-S2)*N1+i1-S1)*N0+i0-S0] =v; }
    HDV inline void add(int i0, int i1, int i2, int i3, int i4, int i5, type v)  {  A[(((((i5-S5)*N4+i4-S4)*N3+i3-S3)*N2+i2-S2)*N1+i1-S1)*N0+i0-S0]+=v; }
    HDV inline void mul(int i0, int i1, int i2, int i3, int i4, int i5, type v)  {  A[(((((i5-S5)*N4+i4-S4)*N3+i3-S3)*N2+i2-S2)*N1+i1-S1)*N0+i0-S0]*=v; }
    HDV inline type&  operator()(int i0, int i1, int i2, int i3, int i4, int i5) {  return A[(((((i5-S5)*N4+i4-S4)*N3+i3-S3)*N2+i2-S2)*N1+i1-S1)*N0+i0-S0]; }
    HDV inline void set(type  v)  {  for(int i=0;i<N5*N4*N3*N2*N1*N0;i++) A[i] = v; }
    HDV inline void set(type *B)  {  for(int i=0;i<N5*N4*N3*N2*N1*N0;i++) A[i] = B[i]; }
    HDV inline int size() {  return N0*N1*N2*N3*N4*N5; }
    HDV inline int length(int d) {  if(d==0) return N0; if(d==1) return N1; if(d==2) return N2; if(d==3) return N3;
                                    if(d==4) return N4; if(d==5) return N5; }
    HDV inline type& at(int p) {  return A[p]; }
    HDV inline index6D indices(int p)
    {
        int Na = N5*N4*N3*N2*N1*N0;
        int ia = p;
        index6D ndx;
        ia %= Na; Na /= N5; ndx.i5 = (int) (ia/Na) + S5;
        ia %= Na; Na /= N4; ndx.i4 = (int) (ia/Na) + S4;
        ia %= Na; Na /= N3; ndx.i3 = (int) (ia/Na) + S3;
        ia %= Na; Na /= N2; ndx.i2 = (int) (ia/Na) + S2;
        ia %= Na; Na /= N1; ndx.i1 = (int) (ia/Na) + S1;
        ia %= Na; Na /= N0; ndx.i0 = (int) (ia/Na) + S0;
        return ndx;
    }

    HDV kFArray5D<type>& subMlast(int i5);  // make 5D array from last  index i5
    HDV kFArray5D<type>& subMfrst(int i0);  // make 5D array from first index i0

    HDV void make(int n0, int n1, int n2, int n3, int n4, int n5)
    {   Z0=0;
        S0= 0; N0= n0; E0= n0-1;
        S1= 0; N1= n1; E1= n1-1;
        S2= 0; N2= n2; E2= n2-1;
        S3= 0; N3= n3; E3= n3-1;
        S4= 0; N4= n4; E4= n4-1;
        S5= 0; N5= n5; E5= n5-1;
        if(A) free(A);
        A = (type *)malloc(N5*N4*N3*N2*N1*N0*sizeof(type));
    }
    HDV void make(int L0, int U0,
                                  int L1, int U1,
                                  int L2, int U2,
                                  int L3, int U3,
                                  int L4, int U4,
                                  int L5, int U5)
    {   Z0= 0;
        S0=L0; N0=U0-L0+1; E0=U0;
        S1=L1; N1=U1-L1+1; E1=U1;
        S2=L2; N2=U2-L2+1; E2=U2;
        S3=L3; N3=U3-L3+1; E3=U3;
        S4=L4; N4=U4-L4+1; E4=U4;
        S5=L5; N5=U5-L5+1; E5=U5;
        if(hA) free(hA); hA = (type *)malloc(N5*N4*N3*N2*N1*N0*sizeof(type));
        dA = NULL; A = hA;
    }

    HST void write(FILE* f)
    {
        fwrite(&Z0,sizeof(int),1,f);
        fwrite(&S0,sizeof(int),1,f);
        fwrite(&N0,sizeof(int),1,f);
        fwrite(&S1,sizeof(int),1,f);
        fwrite(&N1,sizeof(int),1,f);
        fwrite(&S2,sizeof(int),1,f);
        fwrite(&N2,sizeof(int),1,f);
        fwrite(&S3,sizeof(int),1,f);
        fwrite(&N3,sizeof(int),1,f);
        fwrite(&S4,sizeof(int),1,f);
        fwrite(&N4,sizeof(int),1,f);
        fwrite(&S5,sizeof(int),1,f);
        fwrite(&N5,sizeof(int),1,f);
        fwrite(hA,sizeof(type),N0*N1*N2*N3*N4*N5,f);
    }
    HST void read (FILE* f)
    {   if(hA) free(hA);
        fread (&Z0,sizeof(int),1,f);
        fread (&S0,sizeof(int),1,f);
        fread (&N0,sizeof(int),1,f);
        fread (&S1,sizeof(int),1,f);
        fread (&N1,sizeof(int),1,f);
        fread (&S2,sizeof(int),1,f);
        fread (&N2,sizeof(int),1,f);
        fread (&S3,sizeof(int),1,f);
        fread (&N3,sizeof(int),1,f);
        fread (&S4,sizeof(int),1,f);
        fread (&N4,sizeof(int),1,f);
        fread (&S5,sizeof(int),1,f);
        fread (&N5,sizeof(int),1,f);
        hA = (type *)malloc(N0*N1*N2*N3*N4*N5*sizeof(type));
        fread(hA,sizeof(type),N0*N1*N2*N3*N4*N5,f);
        A = hA;
        E0=S0+N0-1;
        E1=S1+N1-1;
        E2=S2+N2-1;
        E3=S3+N3-1;
        E4=S4+N4-1;
        E5=S5+N5-1;
    }

};

template<class type> class kFArray7D
{ // col-major storage: first index varies first
public:
    int Z0;
    int N0,N1,N2,N3,N4,N5,N6;
    int S0,S1,S2,S3,S4,S5,S6;
    int E0,E1,E2,E3,E4,E5,E6;
    type* hA;
    type* dA;
    type* A;
    HDV ~kFArray7D() { if(hA) free(hA); }
    HDV  kFArray7D():Z0(0),
                     S0(0),S1(0),S2(0),S3(0),S4(0),S5(0),S6(0),
                     N0(0),N1(0),N2(0),N3(0),N4(0),N5(0),N6(0),
                     E0(0),E1(0),E2(0),E3(0),E4(0),E5(0),E6(0)
    {  hA = NULL; dA = NULL; A = NULL; }
    HDV  kFArray7D(int n0, int n1, int n2, int n3, int n4, int n5, int n6): Z0(0),
                                   S0(0)   ,S1(0)   ,S2(0)   ,S3(0)   ,S4(0)   ,S5(0)   ,S6(0),
                                   N0(n0  ),N1(n1  ),N2(n2  ),N3(n3  ),N4(n4  ),N5(n5  ),N6(n6  ),
                                   E0(n0-1),E1(n1-1),E2(n2-1),E3(n3-1),E4(n4-1),E5(n5-1),E6(n6-1)
    {  hA = (type *)malloc(N0*N1*N2*N3*N4*N5*N6*sizeof(type)); dA = NULL; A = hA; }
    HDV  kFArray7D(int L0, int U0,
                   int L1, int U1,
                   int L2, int U2,
                   int L3, int U3,
                   int L4, int U4,
                   int L5, int U5,
                   int L6, int U6):Z0(0),
                                   S0(L0)     ,S1(L1)     ,S2(L2)     ,S3(L3)     ,S4(L4)     ,S5(L5)     ,S6(L6),
                                   N0(U0-L0+1),N1(U1-L1+1),N2(U2-L2+1),N3(U3-L3+1),N4(U4-L4+1),N5(U5-L5+1),N6(U6-L6+1),
                                   E0(U0)     ,E1(U1)     ,E2(U2)     ,E3(U3)     ,E4(U4)     ,E5(U5)     ,E6(U6)
    {  hA = (type *)malloc(N0*N1*N2*N3*N4*N5*N6*sizeof(type)); dA = NULL; A = hA; }

    HST kFArray7D<type>* clone()
    {   // copy structure without new allocation
        kFArray7D<type>* hM = new kFArray7D();
        hM->Z0 = Z0;
        hM->S0 = S0;    hM->E0 = E0;    hM->N0 = N0;
        hM->S1 = S1;    hM->E1 = E1;    hM->N1 = N1;
        hM->S2 = S2;    hM->E2 = E2;    hM->N2 = N2;
        hM->S3 = S3;    hM->E3 = E3;    hM->N3 = N3;
        hM->S4 = S4;    hM->E4 = E4;    hM->N4 = N4;
        hM->S5 = S5;    hM->E5 = E5;    hM->N5 = N5;
        hM->S6 = S6;    hM->E6 = E6;    hM->N6 = N6;
        hM->hA = hA;    hM->A  =  A;
        return hM;
    }
#if NVCC
    HST kFArray7D<type>* makeDeviceCopy()
    {
        int Nelmts = N0*N1*N2*N3*N4*N5*N6;
        type* saved_hA = hA;
        if(dA) cudaFree(dA); cudaMalloc((void **) &dA,Nelmts*sizeof(type));
        cudaMemcpy(dA, hA, Nelmts*sizeof(type), cudaMemcpyHostToDevice);
        hA = NULL; A = dA;
        kFArray7D<type>* dM;
        cudaMalloc((void **) &dM,sizeof(*this));
        cudaMemcpy(dM, this, sizeof(*this), cudaMemcpyHostToDevice);
        hA = saved_hA; A = hA;
        return dM;
    }
    HST kFArray7D<type>* makeDeviceOnly()
    {
        int Nelmts = N0*N1*N2*N3*N4*N5*N6;
        type* saved_hA = hA;
        if(dA) cudaFree(dA); cudaMalloc((void **) &dA,Nelmts*sizeof(type));
        // no copying of Array, only pointer info
        hA = NULL; A = dA;
        kFArray7D<type>* dM;
        cudaMalloc((void **) &dM,sizeof(*this));
        cudaMemcpy(dM, this, sizeof(*this), cudaMemcpyHostToDevice);
        hA = saved_hA; A = hA;
        return dM;
    }
    HST kFArray7D<type>* makeDeviceCopy(bool cpyArray)
    {
        if(cpyArray) return makeDeviceCopy();
        else         return makeDeviceOnly();
    }
    HST void freeDeviceCopy()  {  A=hA; if(dA) cudaFree(dA); dA = NULL; }
    HST void freeDeviceCopy(kFArray7D<type>* dM)  {  A=hA; if(dA) cudaFree(dA); dA = NULL; cudaFree(dM); }
    HST void copyDataHostToDevice() { cudaMemcpy(dA, hA, N0*N1*N2*N3*N4*N5*N6*sizeof(type), cudaMemcpyHostToDevice); }
    HST void copyDataDeviceToHost() { cudaMemcpy(hA, dA, N0*N1*N2*N3*N4*N5*N6*sizeof(type), cudaMemcpyDeviceToHost); }
#endif
    HDV inline type get(int i0, int i1, int i2, int i3, int i4, int i5, int i6)   {  return A[((((((i6-S6)*N5+i5-S5)*N4+i4-S4)*N3+i3-S3)*N2+i2-S2)*N1+i1-S1)*N0+i0-S0]; }
    HDV inline void set(int i0, int i1, int i2, int i3, int i4, int i5, int i6, type v)  {  A[((((((i6-S6)*N5+i5-S5)*N4+i4-S4)*N3+i3-S3)*N2+i2-S2)*N1+i1-S1)*N0+i0-S0] =v; }
    HDV inline void add(int i0, int i1, int i2, int i3, int i4, int i5, int i6, type v)  {  A[((((((i6-S6)*N5+i5-S5)*N4+i4-S4)*N3+i3-S3)*N2+i2-S2)*N1+i1-S1)*N0+i0-S0]+=v; }
    HDV inline void mul(int i0, int i1, int i2, int i3, int i4, int i5, int i6, type v)  {  A[((((((i6-S6)*N5+i5-S5)*N4+i4-S4)*N3+i3-S3)*N2+i2-S2)*N1+i1-S1)*N0+i0-S0]*=v; }
    HDV inline type& operator()(int i0, int i1, int i2, int i3, int i4, int i5, int i6) {  return A[((((((i6-S6)*N5+i5-S5)*N4+i4-S4)*N3+i3-S3)*N2+i2-S2)*N1+i1-S1)*N0+i0-S0]; }
    HDV inline void set(type  v)  {  for(int i=0;i<N6*N5*N4*N3*N2*N1*N0;i++) A[i] = v; }
    HDV inline void set(type *B)  {  for(int i=0;i<N6*N5*N4*N3*N2*N1*N0;i++) A[i] = B[i]; }
    HDV inline int size() {  return N0*N1*N2*N3*N4*N5*N6; }
    HDV inline int length(int d) {  if(d==0) return N0; if(d==1) return N1; if(d==2) return N2; if(d==3) return N3;
                                    if(d==4) return N4; if(d==5) return N5; if(d==6) return N6; }
    HDV inline type& at(int p) {  return A[p]; }
    HDV inline index7D indices(int p)
    {
        int Na = N6*N5*N4*N3*N2*N1*N0;
        int ia = p;
        index7D ndx;
        ia %= Na; Na /= N6; ndx.i6 = (int) (ia/Na) + S6;
        ia %= Na; Na /= N5; ndx.i5 = (int) (ia/Na) + S5;
        ia %= Na; Na /= N4; ndx.i4 = (int) (ia/Na) + S4;
        ia %= Na; Na /= N3; ndx.i3 = (int) (ia/Na) + S3;
        ia %= Na; Na /= N2; ndx.i2 = (int) (ia/Na) + S2;
        ia %= Na; Na /= N1; ndx.i1 = (int) (ia/Na) + S1;
        ia %= Na; Na /= N0; ndx.i0 = (int) (ia/Na) + S0;
        return ndx;
    }

    HDV kFArray6D<type>& subMlast(int i6);  // make 6D array from last  index i6
    HDV kFArray6D<type>& subMfrst(int i0);  // make 6D array from first index i0

    HDV void make(int n0, int n1, int n2, int n3, int n4, int n5, int n6)
    {   Z0= 0;
        S0= 0; N0= n0; E0= n0-1;
        S1= 0; N1= n1; E1= n1-1;
        S2= 0; N2= n2; E2= n2-1;
        S3= 0; N3= n3; E3= n3-1;
        S4= 0; N4= n4; E4= n4-1;
        S5= 0; N5= n5; E5= n5-1;
        S6= 0; N6= n6; E6= n6-1;
        if(hA) free(hA); hA = (type *)malloc(N6*N5*N4*N3*N2*N1*N0*sizeof(type));
        dA = NULL; A = hA;
    }
    HDV void make(int L0, int U0,
                                  int L1, int U1,
                                  int L2, int U2,
                                  int L3, int U3,
                                  int L4, int U4,
                                  int L5, int U5,
                                  int L6, int U6)
    {   Z0= 0;
        S0=L0; N0=U0-L0+1; E0=U0;
        S1=L1; N1=U1-L1+1; E1=U1;
        S2=L2; N2=U2-L2+1; E2=U2;
        S3=L3; N3=U3-L3+1; E3=U3;
        S4=L4; N4=U4-L4+1; E4=U4;
        S5=L5; N5=U5-L5+1; E5=U5;
        S6=L6; N6=U6-L6+1; E6=U6;
        if(hA) free(hA); hA = (type *)malloc(N6*N5*N4*N3*N2*N1*N0*sizeof(type));
        dA = NULL; A = hA;
    }

    HST void write(FILE* f)
    {
        fwrite(&Z0,sizeof(int),1,f);
        fwrite(&S0,sizeof(int),1,f);
        fwrite(&N0,sizeof(int),1,f);
        fwrite(&S1,sizeof(int),1,f);
        fwrite(&N1,sizeof(int),1,f);
        fwrite(&S2,sizeof(int),1,f);
        fwrite(&N2,sizeof(int),1,f);
        fwrite(&S3,sizeof(int),1,f);
        fwrite(&N3,sizeof(int),1,f);
        fwrite(&S4,sizeof(int),1,f);
        fwrite(&N4,sizeof(int),1,f);
        fwrite(&S5,sizeof(int),1,f);
        fwrite(&N5,sizeof(int),1,f);
        fwrite(&S6,sizeof(int),1,f);
        fwrite(&N6,sizeof(int),1,f);
        fwrite(hA,sizeof(type),N0*N1*N2*N3*N4*N5*N6,f);
    }
    HST void read (FILE* f)
    {   if(hA) free(hA);
        fread (&Z0,sizeof(int),1,f);
        fread (&S0,sizeof(int),1,f);
        fread (&N0,sizeof(int),1,f);
        fread (&S1,sizeof(int),1,f);
        fread (&N1,sizeof(int),1,f);
        fread (&S2,sizeof(int),1,f);
        fread (&N2,sizeof(int),1,f);
        fread (&S3,sizeof(int),1,f);
        fread (&N3,sizeof(int),1,f);
        fread (&S4,sizeof(int),1,f);
        fread (&N4,sizeof(int),1,f);
        fread (&S5,sizeof(int),1,f);
        fread (&N5,sizeof(int),1,f);
        fread (&S6,sizeof(int),1,f);
        fread (&N6,sizeof(int),1,f);
        hA = (type *)malloc(N0*N1*N2*N3*N4*N5*N6*sizeof(type));
        fread(hA,sizeof(type),N0*N1*N2*N3*N4*N5*N6,f);
        A= hA;
        E0=S0+N0-1;
        E1=S1+N1-1;
        E2=S2+N2-1;
        E3=S3+N3-1;
        E4=S4+N4-1;
        E5=S5+N5-1;
        E6=S6+N6-1;
    }

};


/*
 * C-like arrays: Row-major format
 */

template<class type> class kCArray1D
{ // row-major storage: last index varies first
public:
    int Z0;
    int N0;
    int S0;
    int E0;
    type* hA;
    type* dA;
    type* A;
    HDV ~kCArray1D() { if(hA) free(hA); }
    HDV  kCArray1D(): Z0(0),S0(0),N0(0),E0(0)
    {  hA = NULL; dA = NULL; A = NULL; }
    HDV  kCArray1D(int n0): Z0(0),S0(0),N0(n0),E0(n0-1)
    {  hA = (type *)malloc(N0*sizeof(type)); dA = NULL; A = hA; }
    HDV  kCArray1D(int L0, int U0): Z0(0),S0(L0),N0(U0-L0+1),E0(U0)
    {  hA = (type *)malloc(N0*sizeof(type)); dA = NULL; A = hA; }

    HST kCArray1D<type>* clone()
    {   // copy structure without new allocation
        kCArray1D<type>* hM = new kCArray1D();
        hM->Z0 = Z0;
        hM->S0 = S0;    hM->E0 = E0;    hM->N0 = N0;
        hM->hA = hA;    hM->A  =  A;
        return hM;
    }
#if NVCC
    HST kCArray1D<type>* makeDeviceCopy()
    {
        int Nelmts = N0;
        type* saved_hA = hA;
        if(dA) cudaFree(dA); cudaMalloc((void **) &dA,Nelmts*sizeof(type));
        cudaMemcpy(dA, hA, Nelmts*sizeof(type), cudaMemcpyHostToDevice);
        hA = NULL; A = dA;
        kCArray1D<type>* dM;
        cudaMalloc((void **) &dM,sizeof(*this));
        cudaMemcpy(dM, this, sizeof(*this), cudaMemcpyHostToDevice);
        hA = saved_hA; A = hA;
        return dM;
    }
    HST kCArray1D<type>* makeDeviceOnly()
    {
        int Nelmts = N0;
        type* saved_hA = hA;
        if(dA) cudaFree(dA); cudaMalloc((void **) &dA,Nelmts*sizeof(type));
        // no copying of Array, only pointer info
        hA = NULL; A = dA;
        kCArray1D<type>* dM;
        cudaMalloc((void **) &dM,sizeof(*this));
        cudaMemcpy(dM, this, sizeof(*this), cudaMemcpyHostToDevice);
        hA = saved_hA; A = hA;
        return dM;
    }
    HST kCArray1D<type>* makeDeviceCopy(bool cpyArray)
    {
        if(cpyArray) return makeDeviceCopy();
        else         return makeDeviceOnly();
    }
    HST void freeDeviceCopy()  {  A=hA; if(dA) cudaFree(dA); dA = NULL; }
    HST void freeDeviceCopy(kCArray1D<type>* dM)  {  A=hA; if(dA) cudaFree(dA); dA = NULL; cudaFree(dM); }
    HST void copyDataHostToDevice() { cudaMemcpy(dA, hA, N0*sizeof(type), cudaMemcpyHostToDevice); }
    HST void copyDataDeviceToHost() { cudaMemcpy(hA, dA, N0*sizeof(type), cudaMemcpyDeviceToHost); }
#endif
    HDV  inline type get(int i0)   {  return A[i0-S0]; }
    HDV  inline void set(int i0, type v)  {  A[i0-S0] =v; }
    HDV  inline void add(int i0, type v)  {  A[i0-S0]+=v; }
    HDV  inline void mul(int i0, type v)  {  A[i0-S0]*=v; }
    HDV  inline type&  operator()(int i0) {  return A[i0-S0]; }
    HDV  inline void set(type  v)  {  for(int i=0;i<N0;i++) A[i] = v; }
    HDV  inline void set(type *B)  {  for(int i=0;i<N0;i++) A[i] = B[i]; }
    HDV inline int size()  {  return N0; }
    HDV inline int length() { return N0; }
    HDV inline int length(int d) { return N0; }
    HDV inline type& at(int p) {  return A[p]; }
    HDV inline index1D indices(int p)
    {
        int Na = N0;
        int ia = p;
        index1D ndx;
        ia %= Na; Na /= N0; ndx.i0 = (int) (ia/Na) + S0;
        return ndx;
    }

    HDV void make(int n0)
    {   Z0= 0; S0= 0; N0= n0; E0= n0-1;
        if(hA) free(hA); hA = (type *)malloc(N0*sizeof(type));
        dA = NULL; A = hA;
    }
    HDV void make(int L0, int U0)
    {   Z0= 0; S0=L0; N0=U0-L0+1; E0=U0;
        if(hA) free(hA); hA = (type *)malloc(N0*sizeof(type));
        dA = NULL; A = hA;
    }
    HST void write(FILE* f)
    {
        fwrite(&Z0,sizeof(int),1,f);
        fwrite(&S0,sizeof(int),1,f);
        fwrite(&N0,sizeof(int),1,f);
        fwrite(hA,sizeof(type),N0,f);
    }
    HST void read (FILE* f)
    {   if(hA) free(hA);
        fread (&Z0,sizeof(int),1,f);
        fread (&S0,sizeof(int),1,f);
        fread (&N0,sizeof(int),1,f);
        hA = (type *)malloc(N0*sizeof(type));
        fread(hA,sizeof(type),N0,f);
        E0=S0+N0-1;
    }

};

template<class type> class kCArray2D
{ // row-major storage: last index varies first
public:
    int Z0;
    int N0,N1;  // array size in each direction
    int S0,S1;  // index shift in each direction
    int E0,E1;  // maximum-index in each direction
    type* hA;
    type* dA;
    type* A;
    HDV ~kCArray2D() { if(hA) free(hA); }
    HDV  kCArray2D(): Z0(0),
                      S0(0),S1(0),
                      N0(0),N1(0),
                      E0(0),E1(0)
    {  hA = NULL; dA = NULL; A = NULL; }
    HDV  kCArray2D(int n0, int n1):Z0(0),
                                   S0(0),   S1(0),
                                   N0(n0),  N1(n1),
                                   E0(n0-1),E1(n1-1)
    {  hA = (type *)malloc(N0*N1*sizeof(type)); dA = NULL; A = hA; }
    HDV  kCArray2D(int L0, int U0,
                   int L1, int U1):Z0(0),
                                   S0(L0),     S1(L1),
                                   N0(U0-L0+1),N1(U1-L1+1),
                                   E0(U0),     E1(U1)
    {  hA = (type *)malloc(N0*N1*sizeof(type)); dA = NULL; A = hA; }

    HST kCArray2D<type>* clone()
    {   // copy structure without new allocation
        kCArray2D<type>* hM = new kCArray2D();
        hM->Z0 = Z0;
        hM->S0 = S0;    hM->E0 = E0;    hM->N0 = N0;
        hM->S1 = S1;    hM->E1 = E1;    hM->N1 = N1;
        hM->hA = hA;    hM->A  =  A;
        return hM;
    }
#if NVCC
    HST kCArray2D<type>* makeDeviceCopy()
    {
        int Nelmts = N0*N1;
        type* saved_hA = hA;
        if(dA) cudaFree(dA); cudaMalloc((void **) &dA,Nelmts*sizeof(type));
        cudaMemcpy(dA, hA, Nelmts*sizeof(type), cudaMemcpyHostToDevice);
        hA = NULL; A = dA;
        kCArray2D<type>* dM;
        cudaMalloc((void **) &dM,sizeof(*this));
        cudaMemcpy(dM, this, sizeof(*this), cudaMemcpyHostToDevice);
        hA = saved_hA; A = hA;
        return dM;
    }
    HST kCArray2D<type>* makeDeviceOnly()
    {
        int Nelmts = N0*N1;
        type* saved_hA = hA;
        if(dA) cudaFree(dA); cudaMalloc((void **) &dA,Nelmts*sizeof(type));
        // no copying of Array, only pointer info
        hA = NULL; A = dA;
        kCArray2D<type>* dM;
        cudaMalloc((void **) &dM,sizeof(*this));
        cudaMemcpy(dM, this, sizeof(*this), cudaMemcpyHostToDevice);
        hA = saved_hA; A = hA;
        return dM;
    }
    HST kCArray2D<type>* makeDeviceCopy(bool cpyArray)
    {
        if(cpyArray) return makeDeviceCopy();
        else         return makeDeviceOnly();
    }
    HST void freeDeviceCopy()  {  A=hA; if(dA) cudaFree(dA); dA = NULL; }
    HST void freeDeviceCopy(kCArray2D<type>* dM)  {  A=hA; if(dA) cudaFree(dA); dA = NULL; cudaFree(dM); }
    HST void copyDataHostToDevice() { cudaMemcpy(dA, hA, N0*N1*sizeof(type), cudaMemcpyHostToDevice); }
    HST void copyDataDeviceToHost() { cudaMemcpy(hA, dA, N0*N1*sizeof(type), cudaMemcpyDeviceToHost); }
#endif
    HDV inline type get(int i0, int i1)   {  return A[(i0-S0)*N1+i1-S1]; }
    HDV inline void set(int i0, int i1, type v)  {  A[(i0-S0)*N1+i1-S1] =v; }
    HDV inline void add(int i0, int i1, type v)  {  A[(i0-S0)*N1+i1-S1]+=v; }
    HDV inline void mul(int i0, int i1, type v)  {  A[(i0-S0)*N1+i1-S1]*=v; }
    HDV inline type&  operator()(int i0, int i1) {  return A[(i0-S0)*N1+i1-S1]; }
    HDV inline void set(type  v)  {  for(int i=0;i<N1*N0;i++) A[i] = v; }
    HDV inline void set(type *B)  {  for(int i=0;i<N1*N0;i++) A[i] = B[i]; }
    HDV inline int size() {  return N0*N1; }
    HDV inline int length(int d) {  if(d==0) return N0; if(d==1) return N1; }
    HDV inline type& at(int p) {  return A[p]; }
    HDV inline index2D indices(int p)
    {
        int Na = N1*N0;
        int ia = p;
        index2D ndx;
        ia %= Na; Na /= N0; ndx.i0 = (int) (ia/Na) + S0;
        ia %= Na; Na /= N1; ndx.i1 = (int) (ia/Na) + S1;
        return ndx;
    }

    HDV kCArray1D<type>& subMlast(int i1);  // make 1D array from last  index i1
    HDV kCArray1D<type>& subMfrst(int i0);  // make 1D array from first index i0

    HDV void make(int n0, int n1)
    {   Z0= 0;
        S0= 0; N0= n0; E0= n0-1;
        S1= 0; N1= n1; E1= n1-1;
        if(hA) free(hA); hA = (type *)malloc(N1*N0*sizeof(type));
        dA = NULL; A = hA;
    }
    HDV void make(int L0, int U0,
                                  int L1, int U1)
    {   Z0= 0;
        S0=L0; N0=U0-L0+1; E0=U0;
        S1=L1; N1=U1-L1+1; E1=U1;
        if(hA) free(hA); hA = (type *)malloc(N1*N0*sizeof(type));
        dA = NULL; A = hA;
    }
    HST void write(FILE* f)
    {
        fwrite(&Z0,sizeof(int),1,f);
        fwrite(&S0,sizeof(int),1,f);
        fwrite(&N0,sizeof(int),1,f);
        fwrite(&S1,sizeof(int),1,f);
        fwrite(&N1,sizeof(int),1,f);
        fwrite(hA,sizeof(type),N0*N1,f);
    }
    HST void read (FILE* f)
    {   if(hA) free(hA);
        fread (&Z0,sizeof(int),1,f);
        fread (&S0,sizeof(int),1,f);
        fread (&N0,sizeof(int),1,f);
        fread (&S1,sizeof(int),1,f);
        fread (&N1,sizeof(int),1,f);
        hA = (type *)malloc(N0*N1*sizeof(type));
        fread(hA,sizeof(type),N0*N1,f);
        E0=S0+N0-1;
        E1=S1+N1-1;
    }
};


template<class type> class kCArray3D
{ // row-major storage: last index varies first
public:
    int Z0;
    int N0,N1,N2;
    int S0,S1,S2;
    int E0,E1,E2;
    type* hA;
    type* dA;
    type* A;
    HDV ~kCArray3D() { if(hA) free(hA); }
    HDV  kCArray3D():Z0(0),
                     S0(0),S1(0),S2(0),
                     N0(0),N1(0),N2(0),
                     E0(0),E1(0),E2(0)
    {  hA = NULL; dA = NULL; A = NULL; }
    HDV  kCArray3D(int n0, int n1, int n2): Z0(0),
                                   S0(0),   S1(0),   S2(0),
                                   N0(n0),  N1(n1),  N2(n2),
                                   E0(n0-1),E1(n1-1),E2(n2-1)
    {  hA = (type *)malloc(N0*N1*N2*sizeof(type)); dA = NULL; A = hA; }
    HDV  kCArray3D(int L0, int U0,
                   int L1, int U1,
                   int L2, int U2):Z0(0),
                                   S0(L0),     S1(L1),     S2(L2),
                                   N0(U0-L0+1),N1(U1-L1+1),N2(U2-L2+1),
                                   E0(U0),     E1(U1),     E2(U2)
    {  hA = (type *)malloc(N0*N1*N2*sizeof(type)); dA = NULL; A = hA; }

    HST kCArray3D<type>* clone()
    {   // copy structure without new allocation
        kCArray3D<type>* hM = new kCArray3D();
        hM->Z0 = Z0;
        hM->S0 = S0;    hM->E0 = E0;    hM->N0 = N0;
        hM->S1 = S1;    hM->E1 = E1;    hM->N1 = N1;
        hM->S2 = S2;    hM->E2 = E2;    hM->N2 = N2;
        hM->hA = hA;    hM->A  =  A;
        return hM;
    }
#if NVCC
    HST kCArray3D<type>* makeDeviceCopy()
    {
        int Nelmts = N0*N1*N2;
        type* saved_hA = hA;
        if(dA) cudaFree(dA); cudaMalloc((void **) &dA,Nelmts*sizeof(type));
        cudaMemcpy(dA, hA, Nelmts*sizeof(type), cudaMemcpyHostToDevice);
        hA = NULL; A = dA;
        kCArray3D<type>* dM;
        cudaMalloc((void **) &dM,sizeof(*this));
        cudaMemcpy(dM, this, sizeof(*this), cudaMemcpyHostToDevice);
        hA = saved_hA; A = hA;
        return dM;
    }
    HST kCArray3D<type>* makeDeviceOnly()
    {
        int Nelmts = N0*N1*N2;
        type* saved_hA = hA;
        if(dA) cudaFree(dA); cudaMalloc((void **) &dA,Nelmts*sizeof(type));
        // no copying of Array, only pointer info
        hA = NULL; A = dA;
        kCArray3D<type>* dM;
        cudaMalloc((void **) &dM,sizeof(*this));
        cudaMemcpy(dM, this, sizeof(*this), cudaMemcpyHostToDevice);
        hA = saved_hA; A = hA;
        return dM;
    }
    HST kCArray3D<type>* makeDeviceCopy(bool cpyArray)
    {
        if(cpyArray) return makeDeviceCopy();
        else         return makeDeviceOnly();
    }
    HST void freeDeviceCopy()  {  A=hA; if(dA) cudaFree(dA); dA = NULL; }
    HST void freeDeviceCopy(kCArray3D<type>* dM)  {  A=hA; if(dA) cudaFree(dA); dA = NULL; cudaFree(dM); }
    HST void copyDataHostToDevice() { cudaMemcpy(dA, hA, N0*N1*N2*sizeof(type), cudaMemcpyHostToDevice); }
    HST void copyDataDeviceToHost() { cudaMemcpy(hA, dA, N0*N1*N2*sizeof(type), cudaMemcpyDeviceToHost); }
#endif
    HDV inline type get(int i0, int i1, int i2)   {  return A[((i0-S0)*N1+i1-S1)*N2+i2-S2]; }
    HDV inline void set(int i0, int i1, int i2, type v)  {  A[((i0-S0)*N1+i1-S1)*N2+i2-S2] =v; }
    HDV inline void add(int i0, int i1, int i2, type v)  {  A[((i0-S0)*N1+i1-S1)*N2+i2-S2]+=v; }
    HDV inline void mul(int i0, int i1, int i2, type v)  {  A[((i0-S0)*N1+i1-S1)*N2+i2-S2]*=v; }
    HDV inline type&  operator()(int i0, int i1, int i2) { return A[((i0-S0)*N1+i1-S1)*N2+i2-S2]; }
    HDV inline void set(type  v)  {  for(int i=0;i<N2*N1*N0;i++) A[i] = v; }
    HDV inline void set(type *B)  {  for(int i=0;i<N2*N1*N0;i++) A[i] = B[i]; }
    HDV inline int size() {  return N0*N1*N2; }
    HDV inline int length(int d) {  if(d==0) return N0; if(d==1) return N1; if(d==2) return N2; }
    HDV inline type& at(int p) {  return A[p]; }
    HDV inline index3D indices(int p)
    {
        int Na = N2*N1*N0;
        int ia = p;
        index3D ndx;
        ia %= Na; Na /= N0; ndx.i0 = (int) (ia/Na) + S0;
        ia %= Na; Na /= N1; ndx.i1 = (int) (ia/Na) + S1;
        ia %= Na; Na /= N2; ndx.i2 = (int) (ia/Na) + S2;
        return ndx;
    }

    HDV kCArray2D<type>& subMlast(int i2);  // make 2D array from last  index i2
    HDV kCArray2D<type>& subMfrst(int i0);  // make 2D array from first index i0

    HDV void make(int n0, int n1, int n2)
    {   Z0= 0;
        S0= 0; N0= n0; E0= n0-1;
        S1= 0; N1= n1; E1= n1-1;
        S2= 0; N2= n2; E2= n2-1;
        if(hA) free(hA); hA = (type *)malloc(N2*N1*N0*sizeof(type));
        dA = NULL; A = hA;
    }
    HDV void make(int L0, int U0,
                                  int L1, int U1,
                                  int L2, int U2)
    {   Z0= 0;
        S0=L0; N0=U0-L0+1; E0=U0;
        S1=L1; N1=U1-L1+1; E1=U1;
        S2=L2; N2=U2-L2+1; E2=U2;
        if(hA) free(hA); hA = (type *)malloc(N2*N1*N0*sizeof(type));
        dA = NULL; A = hA;
    }
    HST void write(FILE* f)
    {
        fwrite(&Z0,sizeof(int),1,f);
        fwrite(&S0,sizeof(int),1,f);
        fwrite(&N0,sizeof(int),1,f);
        fwrite(&S1,sizeof(int),1,f);
        fwrite(&N1,sizeof(int),1,f);
        fwrite(&S2,sizeof(int),1,f);
        fwrite(&N2,sizeof(int),1,f);
        fwrite(hA,sizeof(type),N0*N1*N2,f);
    }
    HST void read (FILE* f)
    {   if(hA) free(hA);
        fread (&Z0,sizeof(int),1,f);
        fread (&S0,sizeof(int),1,f);
        fread (&N0,sizeof(int),1,f);
        fread (&S1,sizeof(int),1,f);
        fread (&N1,sizeof(int),1,f);
        fread (&S2,sizeof(int),1,f);
        fread (&N2,sizeof(int),1,f);
        hA = (type *)malloc(N0*N1*N2*sizeof(type));
        fread(hA,sizeof(type),N0*N1*N2,f);
        A = hA;
        E0=S0+N0-1;
        E1=S1+N1-1;
        E2=S2+N2-1;
    }
};

template<class type> class kCArray4D
{ // row-major storage: last index varies first
public:
    int Z0;
    int N0,N1,N2,N3;
    int S0,S1,S2,S3;
    int E0,E1,E2,E3;
    type* hA;
    type* dA;
    type* A;
    HDV  ~kCArray4D() { if(hA) free(hA); }
    HDV   kCArray4D():Z0(0),
                      S0(0),S1(0),S2(0),S3(0),
                      N0(0),N1(0),N2(0),N3(0),
                      E0(0),E1(0),E2(0),E3(0)
    {  hA = NULL; dA = NULL; A = NULL; }
    HDV   kCArray4D(int n0, int n1,
                    int n2, int n3):Z0(0),
                                    S0( 0)  ,S1( 0)  ,S2( 0)  ,S3( 0)  ,
                                    N0(n0  ),N1(n1  ),N2(n2  ),N3(n3  ),
                                    E0(n0-1),E1(n1-1),E2(n2-1),E3(n3-1)
    {  hA = (type *)malloc(N0*N1*N2*N3*sizeof(type)); dA = NULL; A = hA; }
    HDV   kCArray4D(int L0, int U0,
                    int L1, int U1,
                    int L2, int U2,
                    int L3, int U3):Z0(0),
                           S0(L0)     ,S1(L1)     ,S2(L2)     ,S3(L3),
                           N0(U0-L0+1),N1(U1-L1+1),N2(U2-L2+1),N3(U3-L3+1),
                           E0(U0)     ,E1(U1)     ,E2(U2)     ,E3(U3)
    {  hA = (type *)malloc(N0*N1*N2*N3*sizeof(type)); dA = NULL; A = hA; }

    HST kCArray4D<type>* clone()
    {   // copy structure without new allocation
        kCArray4D<type>* hM = new kCArray4D();
        hM->Z0 = Z0;
        hM->S0 = S0;    hM->E0 = E0;    hM->N0 = N0;
        hM->S1 = S1;    hM->E1 = E1;    hM->N1 = N1;
        hM->S2 = S2;    hM->E2 = E2;    hM->N2 = N2;
        hM->S3 = S3;    hM->E3 = E3;    hM->N3 = N3;
        hM->hA = hA;    hM->A  =  A;
        return hM;
    }
#if NVCC
    HST kCArray4D<type>* makeDeviceCopy()
    {
        int Nelmts = N0*N1*N2*N3;
        type* saved_hA = hA;
        if(dA) cudaFree(dA); cudaMalloc((void **) &dA,Nelmts*sizeof(type));
        cudaMemcpy(dA, hA, Nelmts*sizeof(type), cudaMemcpyHostToDevice);
        hA = NULL; A = dA;
        kCArray4D<type>* dM;
        cudaMalloc((void **) &dM,sizeof(*this));
        cudaMemcpy(dM, this, sizeof(*this), cudaMemcpyHostToDevice);
        hA = saved_hA; A = hA;
        return dM;
    }
    HST kCArray4D<type>* makeDeviceOnly()
    {
        int Nelmts = N0*N1*N2*N3;
        type* saved_hA = hA;
        if(dA) cudaFree(dA); cudaMalloc((void **) &dA,Nelmts*sizeof(type));
        // no copying of Array, only pointer info
        hA = NULL; A = dA;
        kCArray4D<type>* dM;
        cudaMalloc((void **) &dM,sizeof(*this));
        cudaMemcpy(dM, this, sizeof(*this), cudaMemcpyHostToDevice);
        hA = saved_hA; A = hA;
        return dM;
    }
    HST kCArray4D<type>* makeDeviceCopy(bool cpyArray)
    {
        if(cpyArray) return makeDeviceCopy();
        else         return makeDeviceOnly();
    }
    HST void freeDeviceCopy()  {  A=hA; if(dA) cudaFree(dA); dA = NULL; }
    HST void freeDeviceCopy(kCArray4D<type>* dM)  {  A=hA; if(dA) cudaFree(dA); dA = NULL; cudaFree(dM); }
    HST void copyDataHostToDevice() { cudaMemcpy(dA, hA, N0*N1*N2*N3*sizeof(type), cudaMemcpyHostToDevice); }
    HST void copyDataDeviceToHost() { cudaMemcpy(hA, dA, N0*N1*N2*N3*sizeof(type), cudaMemcpyDeviceToHost); }
#endif
    HDV inline type get(int i0, int i1, int i2, int i3)   {  return A[(((i0-S0)*N1+i1-S1)*N2+i2-S2)*N3+i3-S3]; }
    HDV inline void set(int i0, int i1, int i2, int i3, type v)  {  A[(((i0-S0)*N1+i1-S1)*N2+i2-S2)*N3+i3-S3] =v; }
    HDV inline void add(int i0, int i1, int i2, int i3, type v)  {  A[(((i0-S0)*N1+i1-S1)*N2+i2-S2)*N3+i3-S3]+=v; }
    HDV inline void mul(int i0, int i1, int i2, int i3, type v)  {  A[(((i0-S0)*N1+i1-S1)*N2+i2-S2)*N3+i3-S3]*=v; }
    HDV inline type&  operator()(int i0, int i1, int i2, int i3) {  return A[(((i0-S0)*N1+i1-S1)*N2+i2-S2)*N3+i3-S3]; }
    HDV inline void set(type  v)  {  for(int i=0;i<N3*N2*N1*N0;i++) A[i] = v; }
    HDV inline void set(type *B)  {  for(int i=0;i<N3*N2*N1*N0;i++) A[i] = B[i]; }
    HDV inline int size() {  return N0*N1*N2*N3; }
    HDV inline int length(int d) {  if(d==0) return N0; if(d==1) return N1; if(d==2) return N2; if(d==3) return N3; }
    HDV inline type& at(int p) {  return A[p]; }
    HDV inline index4D indices(int p)
    {
        int Na = N3*N2*N1*N0;
        int ia = p;
        index4D ndx;
        ia %= Na; Na /= N0; ndx.i0 = (int) (ia/Na) + S0;
        ia %= Na; Na /= N1; ndx.i1 = (int) (ia/Na) + S1;
        ia %= Na; Na /= N2; ndx.i2 = (int) (ia/Na) + S2;
        ia %= Na; Na /= N3; ndx.i3 = (int) (ia/Na) + S3;
        return ndx;
    }

    HDV kCArray3D<type>& subMlast(int i3);  // make 3D array from last  index i3
    HDV kCArray3D<type>& subMfrst(int i0);  // make 3D array from first index i0

    HDV void make(int n0, int n1, int n2, int n3)
    {   Z0= 0;
        S0= 0; N0= n0; E0= n0-1;
        S1= 0; N1= n1; E1= n1-1;
        S2= 0; N2= n2; E2= n2-1;
        S3= 0; N3= n3; E3= n3-1;
        if(hA) free(hA); hA = (type *)malloc(N3*N2*N1*N0*sizeof(type));
        dA = NULL; A = hA;
    }
    HDV void make(int L0, int U0,
                                  int L1, int U1,
                                  int L2, int U2,
                                  int L3, int U3)
    {   Z0= 0;
        S0=L0; N0=U0-L0+1; E0=U0;
        S1=L1; N1=U1-L1+1; E1=U1;
        S2=L2; N2=U2-L2+1; E2=U2;
        S3=L3; N3=U3-L3+1; E3=U3;
        if(hA) free(hA); hA = (type *)malloc(N3*N2*N1*N0*sizeof(type));
        dA = NULL; A = hA;
    }
    HST void write(FILE* f)
    {
        fwrite(&Z0,sizeof(int),1,f);
        fwrite(&S0,sizeof(int),1,f);
        fwrite(&N0,sizeof(int),1,f);
        fwrite(&S1,sizeof(int),1,f);
        fwrite(&N1,sizeof(int),1,f);
        fwrite(&S2,sizeof(int),1,f);
        fwrite(&N2,sizeof(int),1,f);
        fwrite(&S3,sizeof(int),1,f);
        fwrite(&N3,sizeof(int),1,f);
        fwrite(hA,sizeof(type),N0*N1*N2*N3,f);
    }
    HST void read (FILE* f)
    {   if(hA) free(hA);
        fread (&Z0,sizeof(int),1,f);
        fread (&S0,sizeof(int),1,f);
        fread (&N0,sizeof(int),1,f);
        fread (&S1,sizeof(int),1,f);
        fread (&N1,sizeof(int),1,f);
        fread (&S2,sizeof(int),1,f);
        fread (&N2,sizeof(int),1,f);
        fread (&S3,sizeof(int),1,f);
        fread (&N3,sizeof(int),1,f);
        hA = (type *)malloc(N0*N1*N2*N3*sizeof(type));
        fread(hA,sizeof(type),N0*N1*N2*N3,f);
        A = hA;
        E0=S0+N0-1;
        E1=S1+N1-1;
        E2=S2+N2-1;
        E3=S3+N3-1;
    }
};

template<class type> class kCArray5D
{ // row-major storage: last index varies first
public:
    int Z0;
    int N0,N1,N2,N3,N4;
    int S0,S1,S2,S3,S4;
    int E0,E1,E2,E3,E4;
    type* hA;
    type* dA;
    type* A;
    HDV ~kCArray5D() { if(hA) free(hA); }
    HDV  kCArray5D():Z0(0),
                     S0(0),S1(0),S2(0),S3(0),S4(0),
                     N0(0),N1(0),N2(0),N3(0),N4(0),
                     E0(0),E1(0),E2(0),E3(0),E4(0)
    {  hA = NULL; dA = NULL; A = NULL; }
    HDV  kCArray5D(int n0, int n1, int n2, int n3, int n4): Z0(0),
                                   S0(0)   ,S1(0)   ,S2(0)   ,S3(0)   ,S4(0),
                                   N0(n0  ),N1(n1  ),N2(n2  ),N3(n3  ),N4(n4  ),
                                   E0(n0-1),E1(n1-1),E2(n2-1),E3(n3-1),E4(n4-1)
    {  hA = (type *)malloc(N0*N1*N2*N3*N4*sizeof(type)); dA = NULL; A = hA; }
    HDV  kCArray5D(int L0, int U0,
                   int L1, int U1,
                   int L2, int U2,
                   int L3, int U3,
                   int L4, int U4):Z0(0),
                                   S0(L0)     ,S1(L1)     ,S2(L2)     ,S3(L3)     ,S4(L4),
                                   N0(U0-L0+1),N1(U1-L1+1),N2(U2-L2+1),N3(U3-L3+1),N4(U4-L4+1),
                                   E0(U0)     ,E1(U1)     ,E2(U2)     ,E3(U3)     ,E4(U4)
    {  hA = (type *)malloc(N0*N1*N2*N3*N4*sizeof(type)); dA = NULL; A = hA; }

    HST kCArray5D<type>* clone()
    {   // copy structure without new allocation
        kCArray5D<type>* hM = new kCArray5D();
        hM->Z0 = Z0;
        hM->S0 = S0;    hM->E0 = E0;    hM->N0 = N0;
        hM->S1 = S1;    hM->E1 = E1;    hM->N1 = N1;
        hM->S2 = S2;    hM->E2 = E2;    hM->N2 = N2;
        hM->S3 = S3;    hM->E3 = E3;    hM->N3 = N3;
        hM->S4 = S4;    hM->E4 = E4;    hM->N4 = N4;
        hM->hA = hA;    hM->A  =  A;
        return hM;
    }
#if NVCC
    HST kCArray5D<type>* makeDeviceCopy()
    {
        int Nelmts = N0*N1*N2*N3*N4;
        type* saved_hA = hA;
        if(dA) cudaFree(dA); cudaMalloc((void **) &dA,Nelmts*sizeof(type));
        cudaMemcpy(dA, hA, Nelmts*sizeof(type), cudaMemcpyHostToDevice);
        hA = NULL; A = dA;
        kCArray5D<type>* dM;
        cudaMalloc((void **) &dM,sizeof(*this));
        cudaMemcpy(dM, this, sizeof(*this), cudaMemcpyHostToDevice);
        hA = saved_hA; A = hA;
        return dM;
    }
    HST kCArray5D<type>* makeDeviceOnly()
    {
        int Nelmts = N0*N1*N2*N3*N4;
        type* saved_hA = hA;
        if(dA) cudaFree(dA); cudaMalloc((void **) &dA,Nelmts*sizeof(type));
        // no copying of Array, only pointer info
        hA = NULL; A = dA;
        kCArray5D<type>* dM;
        cudaMalloc((void **) &dM,sizeof(*this));
        cudaMemcpy(dM, this, sizeof(*this), cudaMemcpyHostToDevice);
        hA = saved_hA; A = hA;
        return dM;
    }
    HST kCArray5D<type>* makeDeviceCopy(bool cpyArray)
    {
        if(cpyArray) return makeDeviceCopy();
        else         return makeDeviceOnly();
    }
    HST void freeDeviceCopy()  {  A=hA; if(dA) cudaFree(dA); dA = NULL; }
    HST void freeDeviceCopy(kCArray5D<type>* dM)  {  A=hA; if(dA) cudaFree(dA); dA = NULL; cudaFree(dM); }
    HST void copyDataHostToDevice() { cudaMemcpy(dA, hA, N0*N1*N2*N3*N4*sizeof(type), cudaMemcpyHostToDevice); }
    HST void copyDataDeviceToHost() { cudaMemcpy(hA, dA, N0*N1*N2*N3*N4*sizeof(type), cudaMemcpyDeviceToHost); }
#endif
    HDV inline type get(int i0, int i1, int i2, int i3, int i4)   {  return A[((((i0-S0)*N1+i1-S1)*N2+i2-S2)*N3+i3-S3)*N4+i4-S4]; }
    HDV inline void set(int i0, int i1, int i2, int i3, int i4, type v)  {  A[((((i0-S0)*N1+i1-S1)*N2+i2-S2)*N3+i3-S3)*N4+i4-S4] =v; }
    HDV inline void add(int i0, int i1, int i2, int i3, int i4, type v)  {  A[((((i0-S0)*N1+i1-S1)*N2+i2-S2)*N3+i3-S3)*N4+i4-S4]+=v; }
    HDV inline void mul(int i0, int i1, int i2, int i3, int i4, type v)  {  A[((((i0-S0)*N1+i1-S1)*N2+i2-S2)*N3+i3-S3)*N4+i4-S4]*=v; }
    HDV inline type&  operator()(int i0, int i1, int i2, int i3, int i4) {  return A[((((i0-S0)*N1+i1-S1)*N2+i2-S2)*N3+i3-S3)*N4+i4-S4]; }
    HDV inline void set(type  v)  {  for(int i=0;i<N4*N3*N2*N1*N0;i++) A[i] = v; }
    HDV inline void set(type *B)  {  for(int i=0;i<N4*N3*N2*N1*N0;i++) A[i] = B[i]; }
    HDV inline int size() {  return N0*N1*N2*N3*N4; }
    HDV inline int length(int d) {  if(d==0) return N0; if(d==1) return N1; if(d==2) return N2; if(d==3) return N3;
                                    if(d==4) return N4; }
    HDV inline type& at(int p) {  return A[p]; }
    HDV inline index5D indices(int p)
    {
        int Na = N4*N3*N2*N1*N0;
        int ia = p;
        index5D ndx;
        ia %= Na; Na /= N0; ndx.i0 = (int) (ia/Na) + S0;
        ia %= Na; Na /= N1; ndx.i1 = (int) (ia/Na) + S1;
        ia %= Na; Na /= N2; ndx.i2 = (int) (ia/Na) + S2;
        ia %= Na; Na /= N3; ndx.i3 = (int) (ia/Na) + S3;
        ia %= Na; Na /= N4; ndx.i4 = (int) (ia/Na) + S4;
        return ndx;
    }

    HDV kCArray4D<type>& subMlast(int i4);  // make 4D array from last  index i4
    HDV kCArray4D<type>& subMfrst(int i0);  // make 4D array from first index i0

    HDV void make(int n0, int n1, int n2, int n3, int n4)
    {   Z0= 0;
        S0= 0; N0= n0; E0= n0-1;
        S1= 0; N1= n1; E1= n1-1;
        S2= 0; N2= n2; E2= n2-1;
        S3= 0; N3= n3; E3= n3-1;
        S4= 0; N4= n4; E4= n4-1;
        if(hA) free(hA); hA = (type *)malloc(N4*N3*N2*N1*N0*sizeof(type));
        dA = NULL; A = hA;
    }
    HDV void make(int L0, int U0,
                  int L1, int U1,
                  int L2, int U2,
                  int L3, int U3,
                  int L4, int U4)
    {   Z0= 0;
        S0=L0; N0=U0-L0+1; E0=U0;
        S1=L1; N1=U1-L1+1; E1=U1;
        S2=L2; N2=U2-L2+1; E2=U2;
        S3=L3; N3=U3-L3+1; E3=U3;
        S4=L4; N4=U4-L4+1; E4=U4;
        if(hA) free(hA); hA = (type *)malloc(N4*N3*N2*N1*N0*sizeof(type));
        dA = NULL; A = hA;
    }
    HST void write(FILE* f)
    {
        fwrite(&Z0,sizeof(int),1,f);
        fwrite(&S0,sizeof(int),1,f);
        fwrite(&N0,sizeof(int),1,f);
        fwrite(&S1,sizeof(int),1,f);
        fwrite(&N1,sizeof(int),1,f);
        fwrite(&S2,sizeof(int),1,f);
        fwrite(&N2,sizeof(int),1,f);
        fwrite(&S3,sizeof(int),1,f);
        fwrite(&N3,sizeof(int),1,f);
        fwrite(&S4,sizeof(int),1,f);
        fwrite(&N4,sizeof(int),1,f);
        fwrite(hA,sizeof(type),N0*N1*N2*N3*N4,f);
    }
    HST void read (FILE* f)
    {   if(hA) free(hA);
        fread (&Z0,sizeof(int),1,f);
        fread (&S0,sizeof(int),1,f);
        fread (&N0,sizeof(int),1,f);
        fread (&S1,sizeof(int),1,f);
        fread (&N1,sizeof(int),1,f);
        fread (&S2,sizeof(int),1,f);
        fread (&N2,sizeof(int),1,f);
        fread (&S3,sizeof(int),1,f);
        fread (&N3,sizeof(int),1,f);
        fread (&S4,sizeof(int),1,f);
        fread (&N4,sizeof(int),1,f);
        hA = (type *)malloc(N0*N1*N2*N3*N4*sizeof(type));
        fread(hA,sizeof(type),N0*N1*N2*N3*N4,f);
        A = hA;
        E0=S0+N0-1;
        E1=S1+N1-1;
        E2=S2+N2-1;
        E3=S3+N3-1;
        E4=S4+N4-1;
    }

};

template<class type> class kCArray6D
{ // row-major storage: last index varies first
public:
    int Z0;
    int N0,N1,N2,N3,N4,N5;
    int S0,S1,S2,S3,S4,S5;
    int E0,E1,E2,E3,E4,E5;
    type* hA;
    type* dA;
    type* A;
    HDV ~kCArray6D() { if(hA) free(hA); }
    HDV  kCArray6D():Z0(0),
                     S0(0),S1(0),S2(0),S3(0),S4(0),S5(0),
                     N0(0),N1(0),N2(0),N3(0),N4(0),N5(0),
                     E0(0),E1(0),E2(0),E3(0),E4(0),E5(0)
    {  hA = NULL; dA = NULL; A = NULL; }
    HDV  kCArray6D(int n0, int n1, int n2, int n3, int n4, int n5): Z0(0),
                                   S0(0)   ,S1(0)   ,S2(0)   ,S3(0)   ,S4(0)   ,S5(0),
                                   N0(n0  ),N1(n1  ),N2(n2  ),N3(n3  ),N4(n4  ),N5(n5  ),
                                   E0(n0-1),E1(n1-1),E2(n2-1),E3(n3-1),E4(n4-1),E5(n5-1)
    {  hA = (type *)malloc(N0*N1*N2*N3*N4*N5*sizeof(type)); dA = NULL; A = hA; }
    HDV  kCArray6D(int L0, int U0,
                   int L1, int U1,
                   int L2, int U2,
                   int L3, int U3,
                   int L4, int U4,
                   int L5, int U5):Z0(0),
                                   S0(L0)     ,S1(L1)     ,S2(L2)     ,S3(L3)     ,S4(L4)     ,S5(L5),
                                   N0(U0-L0+1),N1(U1-L1+1),N2(U2-L2+1),N3(U3-L3+1),N4(U4-L4+1),N5(U5-L5+1),
                                   E0(U0)     ,E1(U1)     ,E2(U2)     ,E3(U3)     ,E4(U4)     ,E5(U5)
    {  hA = (type *)malloc(N0*N1*N2*N3*N4*N5*sizeof(type)); dA = NULL; A = hA; }

    HST kCArray6D<type>* clone()
    {   // copy structure without new allocation
        kCArray6D<type>* hM = new kCArray6D();
        hM->Z0 = Z0;
        hM->S0 = S0;    hM->E0 = E0;    hM->N0 = N0;
        hM->S1 = S1;    hM->E1 = E1;    hM->N1 = N1;
        hM->S2 = S2;    hM->E2 = E2;    hM->N2 = N2;
        hM->S3 = S3;    hM->E3 = E3;    hM->N3 = N3;
        hM->S4 = S4;    hM->E4 = E4;    hM->N4 = N4;
        hM->S5 = S5;    hM->E5 = E5;    hM->N5 = N5;
        hM->hA = hA;    hM->A  =  A;
        return hM;
    }
#if NVCC
    HST kCArray6D<type>* makeDeviceCopy()
    {
        int Nelmts = N0*N1*N2*N3*N4*N5;
        type* saved_hA = hA;
        if(dA) cudaFree(dA); cudaMalloc((void **) &dA,Nelmts*sizeof(type));
        cudaMemcpy(dA, hA, Nelmts*sizeof(type), cudaMemcpyHostToDevice);
        hA = NULL; A = dA;
        kCArray6D<type>* dM;
        cudaMalloc((void **) &dM,sizeof(*this));
        cudaMemcpy(dM, this, sizeof(*this), cudaMemcpyHostToDevice);
        hA = saved_hA; A = hA;
        return dM;
    }
    HST kCArray6D<type>* makeDeviceOnly()
    {
        int Nelmts = N0*N1*N2*N3*N4*N5;
        type* saved_hA = hA;
        if(dA) cudaFree(dA); cudaMalloc((void **) &dA,Nelmts*sizeof(type));
        // no copying of Array, only pointer info
        hA = NULL; A = dA;
        kCArray6D<type>* dM;
        cudaMalloc((void **) &dM,sizeof(*this));
        cudaMemcpy(dM, this, sizeof(*this), cudaMemcpyHostToDevice);
        hA = saved_hA; A = hA;
        return dM;
    }
    HST kCArray6D<type>* makeDeviceCopy(bool cpyArray)
    {
        if(cpyArray) return makeDeviceCopy();
        else         return makeDeviceOnly();
    }
    HST void freeDeviceCopy()  {  A=hA; if(dA) cudaFree(dA); dA = NULL; }
    HST void freeDeviceCopy(kCArray6D<type>* dM)  {  A=hA; if(dA) cudaFree(dA); dA = NULL; cudaFree(dM); }
    HST void copyDataHostToDevice() { cudaMemcpy(dA, hA, N0*N1*N2*N3*N4*N5*sizeof(type), cudaMemcpyHostToDevice); }
    HST void copyDataDeviceToHost() { cudaMemcpy(hA, dA, N0*N1*N2*N3*N4*N5*sizeof(type), cudaMemcpyDeviceToHost); }
#endif
    HDV inline type get(int i0, int i1, int i2, int i3, int i4, int i5)   {  return A[(((((i0-S0)*N1+i1-S1)*N2+i2-S2)*N3+i3-S3)*N4+i4-S4)*N5+i5-S5]; }
    HDV inline void set(int i0, int i1, int i2, int i3, int i4, int i5, type v)  {  A[(((((i0-S0)*N1+i1-S1)*N2+i2-S2)*N3+i3-S3)*N4+i4-S4)*N5+i5-S5] =v; }
    HDV inline void add(int i0, int i1, int i2, int i3, int i4, int i5, type v)  {  A[(((((i0-S0)*N1+i1-S1)*N2+i2-S2)*N3+i3-S3)*N4+i4-S4)*N5+i5-S5]+=v; }
    HDV inline void mul(int i0, int i1, int i2, int i3, int i4, int i5, type v)  {  A[(((((i0-S0)*N1+i1-S1)*N2+i2-S2)*N3+i3-S3)*N4+i4-S4)*N5+i5-S5]*=v; }
    HDV inline type&  operator()(int i0, int i1, int i2, int i3, int i4, int i5) {  return A[(((((i0-S0)*N1+i1-S1)*N2+i2-S2)*N3+i3-S3)*N4+i4-S4)*N5+i5-S5]; }
    HDV inline void set(type  v)  {  for(int i=0;i<N5*N4*N3*N2*N1*N0;i++) A[i] = v; }
    HDV inline void set(type *B)  {  for(int i=0;i<N5*N4*N3*N2*N1*N0;i++) A[i] = B[i]; }
    HDV inline int size() {  return N0*N1*N2*N3*N4*N5; }
    HDV inline int length(int d) {  if(d==0) return N0; if(d==1) return N1; if(d==2) return N2; if(d==3) return N3;
                                    if(d==4) return N4; if(d==5) return N5; }
    HDV inline type& at(int p) {  return A[p]; }
    HDV inline index6D indices(int p)
    {
        int Na = N5*N4*N3*N2*N1*N0;
        int ia = p;
        index6D ndx;
        ia %= Na; Na /= N0; ndx.i0 = (int) (ia/Na) + S0;
        ia %= Na; Na /= N1; ndx.i1 = (int) (ia/Na) + S1;
        ia %= Na; Na /= N2; ndx.i2 = (int) (ia/Na) + S2;
        ia %= Na; Na /= N3; ndx.i3 = (int) (ia/Na) + S3;
        ia %= Na; Na /= N4; ndx.i4 = (int) (ia/Na) + S4;
        ia %= Na; Na /= N5; ndx.i5 = (int) (ia/Na) + S5;
        return ndx;
    }

    HDV kCArray5D<type>& subMlast(int i5);  // make 5D array from last  index i5
    HDV kCArray5D<type>& subMfrst(int i0);  // make 5D array from first index i0

    HDV void make(int n0, int n1, int n2, int n3, int n4, int n5)
    {   Z0= 0;
        S0= 0; N0= n0; E0= n0-1;
        S1= 0; N1= n1; E1= n1-1;
        S2= 0; N2= n2; E2= n2-1;
        S3= 0; N3= n3; E3= n3-1;
        S4= 0; N4= n4; E4= n4-1;
        S5= 0; N5= n5; E5= n5-1;
        if(hA) free(hA); hA = (type *)malloc(N5*N4*N3*N2*N1*N0*sizeof(type));
        dA = NULL; A = hA;
    }
    HDV void make(int L0, int U0,
                  int L1, int U1,
                  int L2, int U2,
                  int L3, int U3,
                  int L4, int U4,
                  int L5, int U5)
    {   Z0= 0;
        S0=L0; N0=U0-L0+1; E0=U0;
        S1=L1; N1=U1-L1+1; E1=U1;
        S2=L2; N2=U2-L2+1; E2=U2;
        S3=L3; N3=U3-L3+1; E3=U3;
        S4=L4; N4=U4-L4+1; E4=U4;
        S5=L5; N5=U5-L5+1; E5=U5;
        if(hA) free(hA); hA = (type *)malloc(N5*N4*N3*N2*N1*N0*sizeof(type));
        dA = NULL; A = hA;
    }
    HST void write(FILE* f)
    {
        fwrite(&Z0,sizeof(int),1,f);
        fwrite(&S0,sizeof(int),1,f);
        fwrite(&N0,sizeof(int),1,f);
        fwrite(&S1,sizeof(int),1,f);
        fwrite(&N1,sizeof(int),1,f);
        fwrite(&S2,sizeof(int),1,f);
        fwrite(&N2,sizeof(int),1,f);
        fwrite(&S3,sizeof(int),1,f);
        fwrite(&N3,sizeof(int),1,f);
        fwrite(&S4,sizeof(int),1,f);
        fwrite(&N4,sizeof(int),1,f);
        fwrite(&S5,sizeof(int),1,f);
        fwrite(&N5,sizeof(int),1,f);
        fwrite(hA,sizeof(type),N0*N1*N2*N3*N4*N5,f);
    }
    HST void read (FILE* f)
    {   if(hA) free(hA);
        fread (&Z0,sizeof(int),1,f);
        fread (&S0,sizeof(int),1,f);
        fread (&N0,sizeof(int),1,f);
        fread (&S1,sizeof(int),1,f);
        fread (&N1,sizeof(int),1,f);
        fread (&S2,sizeof(int),1,f);
        fread (&N2,sizeof(int),1,f);
        fread (&S3,sizeof(int),1,f);
        fread (&N3,sizeof(int),1,f);
        fread (&S4,sizeof(int),1,f);
        fread (&N4,sizeof(int),1,f);
        fread (&S5,sizeof(int),1,f);
        fread (&N5,sizeof(int),1,f);
        hA = (type *)malloc(N0*N1*N2*N3*N4*N5*sizeof(type));
        fread(hA,sizeof(type),N0*N1*N2*N3*N4*N5,f);
        E0=S0+N0-1;
        E1=S1+N1-1;
        E2=S2+N2-1;
        E3=S3+N3-1;
        E4=S4+N4-1;
        E5=S5+N5-1;
    }

};

template<class type> class kCArray7D
{ // row-major storage: last index varies first
public:
    int Z0;
    int N0,N1,N2,N3,N4,N5,N6;
    int S0,S1,S2,S3,S4,S5,S6;
    int E0,E1,E2,E3,E4,E5,E6;
    type* hA;
    type* dA;
    type* A;
    HDV ~kCArray7D() { if(hA) free(hA); }
    HDV  kCArray7D():Z0(0),
                     S0(0),S1(0),S2(0),S3(0),S4(0),S5(0),S6(0),
                     N0(0),N1(0),N2(0),N3(0),N4(0),N5(0),N6(0),
                     E0(0),E1(0),E2(0),E3(0),E4(0),E5(0),E6(0)
    {  hA = NULL; dA = NULL; A = NULL; }
    HDV  kCArray7D(int n0, int n1, int n2, int n3, int n4, int n5, int n6): Z0(0),
                                   S0(0)   ,S1(0)   ,S2(0)   ,S3(0)   ,S4(0)   ,S5(0)   ,S6(0),
                                   N0(n0  ),N1(n1  ),N2(n2  ),N3(n3  ),N4(n4  ),N5(n5  ),N6(n6  ),
                                   E0(n0-1),E1(n1-1),E2(n2-1),E3(n3-1),E4(n4-1),E5(n5-1),E6(n6-1)
    {  hA = (type *)malloc(N0*N1*N2*N3*N4*N5*N6*sizeof(type)); dA = NULL; A = hA; }
    HDV  kCArray7D(int L0, int U0,
                   int L1, int U1,
                   int L2, int U2,
                   int L3, int U3,
                   int L4, int U4,
                   int L5, int U5,
                   int L6, int U6):Z0(0),
                                   S0(L0)     ,S1(L1)     ,S2(L2)     ,S3(L3)     ,S4(L4)     ,S5(L5)     ,S6(L6),
                                   N0(U0-L0+1),N1(U1-L1+1),N2(U2-L2+1),N3(U3-L3+1),N4(U4-L4+1),N5(U5-L5+1),N6(U6-L6+1),
                                   E0(U0)     ,E1(U1)     ,E2(U2)     ,E3(U3)     ,E4(U4)     ,E5(U5)     ,E6(U6)
    {  hA = (type *)malloc(N0*N1*N2*N3*N4*N5*N6*sizeof(type)); dA = NULL; A = hA; }

    HST kCArray7D<type>* clone()
    {   // copy structure without new allocation
        kCArray7D<type>* hM = new kCArray7D();
        hM->Z0 = Z0;
        hM->S0 = S0;    hM->E0 = E0;    hM->N0 = N0;
        hM->S1 = S1;    hM->E1 = E1;    hM->N1 = N1;
        hM->S2 = S2;    hM->E2 = E2;    hM->N2 = N2;
        hM->S3 = S3;    hM->E3 = E3;    hM->N3 = N3;
        hM->S4 = S4;    hM->E4 = E4;    hM->N4 = N4;
        hM->S5 = S5;    hM->E5 = E5;    hM->N5 = N5;
        hM->S6 = S6;    hM->E6 = E6;    hM->N6 = N6;
        hM->hA = hA;    hM->A  =  A;
        return hM;
    }
#if NVCC
    HST kCArray7D<type>* makeDeviceCopy()
    {
        int Nelmts = N0*N1*N2*N3*N4*N5*N6;
        type* saved_hA = hA;
        if(dA) cudaFree(dA); cudaMalloc((void **) &dA,Nelmts*sizeof(type));
        cudaMemcpy(dA, hA, Nelmts*sizeof(type), cudaMemcpyHostToDevice);
        hA = NULL; A = dA;
        kCArray7D<type>* dM;
        cudaMalloc((void **) &dM,sizeof(*this));
        cudaMemcpy(dM, this, sizeof(*this), cudaMemcpyHostToDevice);
        hA = saved_hA; A = hA;
        return dM;
    }
    HST kCArray7D<type>* makeDeviceOnly()
    {
        int Nelmts = N0*N1*N2*N3*N4*N5*N6;
        type* saved_hA = hA;
        if(dA) cudaFree(dA); cudaMalloc((void **) &dA,Nelmts*sizeof(type));
        // no copying of Array, only pointer info
        hA = NULL; A = dA;
        kCArray7D<type>* dM;
        cudaMalloc((void **) &dM,sizeof(*this));
        cudaMemcpy(dM, this, sizeof(*this), cudaMemcpyHostToDevice);
        hA = saved_hA; A = hA;
        return dM;
    }
    HST kCArray7D<type>* makeDeviceCopy(bool cpyArray)
    {
        if(cpyArray) return makeDeviceCopy();
        else         return makeDeviceOnly();
    }
    HST void freeDeviceCopy()  {  A=hA; if(dA) cudaFree(dA); dA = NULL; }
    HST void freeDeviceCopy(kCArray7D<type>* dM)  {  A=hA; if(dA) cudaFree(dA); dA = NULL; cudaFree(dM); }
    HST void copyDataHostToDevice() { cudaMemcpy(dA, hA, N0*N1*N2*N3*N4*N5*N6*sizeof(type), cudaMemcpyHostToDevice); }
    HST void copyDataDeviceToHost() { cudaMemcpy(hA, dA, N0*N1*N2*N3*N4*N5*N6*sizeof(type), cudaMemcpyDeviceToHost); }
#endif
    HDV inline type get(int i0, int i1, int i2, int i3, int i4, int i5, int i6)   {  return A[((((((i0-S0)*N1+i1-S1)*N2+i2-S2)*N3+i3-S3)*N4+i4-S4)*N5+i5-S5)*N6+i6-S6]; }
    HDV inline void set(int i0, int i1, int i2, int i3, int i4, int i5, int i6, type v)  {  A[((((((i0-S0)*N1+i1-S1)*N2+i2-S2)*N3+i3-S3)*N4+i4-S4)*N5+i5-S5)*N6+i6-S6] =v; }
    HDV inline void add(int i0, int i1, int i2, int i3, int i4, int i5, int i6, type v)  {  A[((((((i0-S0)*N1+i1-S1)*N2+i2-S2)*N3+i3-S3)*N4+i4-S4)*N5+i5-S5)*N6+i6-S6]+=v; }
    HDV inline void mul(int i0, int i1, int i2, int i3, int i4, int i5, int i6, type v)  {  A[((((((i0-S0)*N1+i1-S1)*N2+i2-S2)*N3+i3-S3)*N4+i4-S4)*N5+i5-S5)*N6+i6-S6]*=v; }
    HDV inline type& operator()(int i0, int i1, int i2, int i3, int i4, int i5, int i6) {  return A[((((((i0-S0)*N1+i1-S1)*N2+i2-S2)*N3+i3-S3)*N4+i4-S4)*N5+i5-S5)*N6+i6-S6]; }
    HDV inline void set(type  v)  {  for(int i=0;i<N6*N5*N4*N3*N2*N1*N0;i++) A[i] = v; }
    HDV inline void set(type *B)  {  for(int i=0;i<N6*N5*N4*N3*N2*N1*N0;i++) A[i] = B[i]; }
    HDV inline int size() {  return N0*N1*N2*N3*N4*N5*N6; }
    HDV inline int length(int d) {  if(d==0) return N0; if(d==1) return N1; if(d==2) return N2; if(d==3) return N3;
                                    if(d==4) return N4; if(d==5) return N5; if(d==6) return N6; }
    HDV inline type& at(int p) {  return A[p]; }
    HDV inline index7D indices(int p)
    {
        int Na = N6*N5*N4*N3*N2*N1*N0;
        int ia = p;
        index7D ndx;
        ia %= Na; Na /= N0; ndx.i0 = (int) (ia/Na) + S0;
        ia %= Na; Na /= N1; ndx.i1 = (int) (ia/Na) + S1;
        ia %= Na; Na /= N2; ndx.i2 = (int) (ia/Na) + S2;
        ia %= Na; Na /= N3; ndx.i3 = (int) (ia/Na) + S3;
        ia %= Na; Na /= N4; ndx.i4 = (int) (ia/Na) + S4;
        ia %= Na; Na /= N5; ndx.i5 = (int) (ia/Na) + S5;
        ia %= Na; Na /= N6; ndx.i6 = (int) (ia/Na) + S6;
        return ndx;
    }

    HDV kCArray6D<type>& subMlast(int i6);  // make 6D array from last  index i6
    HDV kCArray6D<type>& subMfrst(int i0);  // make 6D array from first index i0

    HDV void make(int n0, int n1, int n2, int n3, int n4, int n5, int n6)
    {   Z0= 0;
        S0= 0; N0= n0; E0= n0-1;
        S1= 0; N1= n1; E1= n1-1;
        S2= 0; N2= n2; E2= n2-1;
        S3= 0; N3= n3; E3= n3-1;
        S4= 0; N4= n4; E4= n4-1;
        S5= 0; N5= n5; E5= n5-1;
        S6= 0; N6= n6; E6= n6-1;
        if(hA) free(hA); hA = (type *)malloc(N6*N5*N4*N3*N2*N1*N0*sizeof(type));
        dA = NULL; A = hA;
    }
    HDV void make(int L0, int U0,
                  int L1, int U1,
                  int L2, int U2,
                  int L3, int U3,
                  int L4, int U4,
                  int L5, int U5,
                  int L6, int U6)
    {   Z0= 0;
        S0=L0; N0=U0-L0+1; E0=U0;
        S1=L1; N1=U1-L1+1; E1=U1;
        S2=L2; N2=U2-L2+1; E2=U2;
        S3=L3; N3=U3-L3+1; E3=U3;
        S4=L4; N4=U4-L4+1; E4=U4;
        S5=L5; N5=U5-L5+1; E5=U5;
        S6=L6; N6=U6-L6+1; E6=U6;
        if(hA) free(hA); hA = (type *)malloc(N6*N5*N4*N3*N2*N1*N0*sizeof(type));
        dA = NULL; A = hA;
    }
    HST void write(FILE* f)
    {
        fwrite(&Z0,sizeof(int),1,f);
        fwrite(&S0,sizeof(int),1,f);
        fwrite(&N0,sizeof(int),1,f);
        fwrite(&S1,sizeof(int),1,f);
        fwrite(&N1,sizeof(int),1,f);
        fwrite(&S2,sizeof(int),1,f);
        fwrite(&N2,sizeof(int),1,f);
        fwrite(&S3,sizeof(int),1,f);
        fwrite(&N3,sizeof(int),1,f);
        fwrite(&S4,sizeof(int),1,f);
        fwrite(&N4,sizeof(int),1,f);
        fwrite(&S5,sizeof(int),1,f);
        fwrite(&N5,sizeof(int),1,f);
        fwrite(&S6,sizeof(int),1,f);
        fwrite(&N6,sizeof(int),1,f);
        fwrite(hA,sizeof(type),N0*N1*N2*N3*N4*N5*N6,f);
    }
    HST void read (FILE* f)
    {   if(hA) free(hA);
        fread (&Z0,sizeof(int),1,f);
        fread (&S0,sizeof(int),1,f);
        fread (&N0,sizeof(int),1,f);
        fread (&S1,sizeof(int),1,f);
        fread (&N1,sizeof(int),1,f);
        fread (&S2,sizeof(int),1,f);
        fread (&N2,sizeof(int),1,f);
        fread (&S3,sizeof(int),1,f);
        fread (&N3,sizeof(int),1,f);
        fread (&S4,sizeof(int),1,f);
        fread (&N4,sizeof(int),1,f);
        fread (&S5,sizeof(int),1,f);
        fread (&N5,sizeof(int),1,f);
        fread (&S6,sizeof(int),1,f);
        fread (&N6,sizeof(int),1,f);
        hA = (type *)malloc(N0*N1*N2*N3*N4*N5*N6*sizeof(type));
        fread(hA,sizeof(type),N0*N1*N2*N3*N4*N5*N6,f);
        A = hA;
        E0=S0+N0-1;
        E1=S1+N1-1;
        E2=S2+N2-1;
        E3=S3+N3-1;
        E4=S4+N4-1;
        E5=S5+N5-1;
        E6=S6+N6-1;
    }
};

// data with float4

class kF4Array1D
{ // col-major storage: first index varies first
public:
    int Z0;
    int N0;
    int S0;
    int E0;
    float4* hA;
    float4* dA;
    float4* A;
    HDV ~kF4Array1D() { if(hA) free(hA); }
    HDV  kF4Array1D(): Z0(0),S0(0),N0(0),E0(0)
    {  hA = NULL; dA = NULL; A = NULL; }
    HDV  kF4Array1D(int n0): Z0(0), S0(0), N0(n0),E0(n0-1)
    {  hA = (float4 *)malloc(N0*sizeof(float4)); dA = NULL; A = hA; }
    HDV  kF4Array1D(int L0, int U0): Z0(0),S0(L0),N0(U0-L0+1),E0(U0)
    {  hA = (float4 *)malloc(N0*sizeof(float4)); dA = NULL; A = hA; }
#if NVCC
    HST kF4Array1D* makeDeviceCopy()
    {
        int Nelmts = N0;
        float4* saved_hA = hA;
        if(dA) cudaFree(dA); cudaMalloc((void **) &dA,Nelmts*sizeof(float4));
        cudaMemcpy(dA, hA, Nelmts*sizeof(float4), cudaMemcpyHostToDevice);
        hA = NULL; A = dA;
        kF4Array1D* dM;
        cudaMalloc((void **) &dM,sizeof(*this));
        cudaMemcpy(dM, this, sizeof(*this), cudaMemcpyHostToDevice);
        hA = saved_hA; A = hA;
        return dM;
    }
    HST kF4Array1D* makeDeviceOnly()
    {
        int Nelmts = N0;
        float4* saved_hA = hA;
        if(dA) cudaFree(dA); cudaMalloc((void **) &dA,Nelmts*sizeof(float4));
        // no copying of Array, only pointer info
        hA = NULL; A = dA;
        kF4Array1D* dM;
        cudaMalloc((void **) &dM,sizeof(*this));
        cudaMemcpy(dM, this, sizeof(*this), cudaMemcpyHostToDevice);
        hA = saved_hA; A = hA;
        return dM;
    }
    HST kF4Array1D* makeDeviceCopy(bool cpyArray)
    {
        if(cpyArray) return makeDeviceCopy();
        else         return makeDeviceOnly();
    }
    HST void freeDeviceCopy()  {  A=hA; if(dA) cudaFree(dA); dA = NULL; }
    HST void freeDeviceCopy(kF4Array1D* dM)  {  A=hA; if(dA) cudaFree(dA); dA = NULL; cudaFree(dM); }
    HST void copyDataHostToDevice() { cudaMemcpy(dA, hA, N0*sizeof(float4), cudaMemcpyHostToDevice); }
    HST void copyDataDeviceToHost() { cudaMemcpy(hA, dA, N0*sizeof(float4), cudaMemcpyDeviceToHost); }
#endif
    HDV inline float4 get(int i0)   {  return A[i0-S0]; }
    HDV inline void set(int i0, float4 v)  {  A[i0-S0] = v; }
    HDV inline void add(int i0, float4 v)  {  A[i0-S0].x += v.x; A[i0-S0].y += v.y;
                                                              A[i0-S0].z += v.z; A[i0-S0].w += v.w; }
    HDV inline void mul(int i0, float4 v)  {  A[i0-S0].x *= v.x; A[i0-S0].y *= v.y;
                                                              A[i0-S0].z *= v.z; A[i0-S0].w *= v.w; }
    HDV inline float4&  operator()(int i0) {  return A[i0-S0]; }
    HDV inline void set(float4  v)  {  for(int i=0;i<N0;i++) A[i] = v; }
    HDV inline void set(float4* B)  {  for(int i=0;i<N0;i++) A[i] = B[i]; }
    HDV inline int size() {  return N0; }
    HDV inline float4& at(int p) {  return A[p]; }

    HDV void make(int n0)
    {   Z0= 0; S0= 0; N0= n0; E0= n0-1;
        if(hA) free(hA); hA = (float4 *)malloc(N0*sizeof(float4));
        dA = NULL; A = hA;
    }
    HDV void make(int L0, int U0)
    {   Z0= 0; S0=L0; N0=U0-L0+1; E0=U0;
        if(hA) free(hA); hA = (float4 *)malloc(N0*sizeof(float4));
        dA = NULL; A = hA;
    }

    HST void write(FILE* f)
    {
        fwrite(&Z0,sizeof(int),1,f);
        fwrite(&S0,sizeof(int),1,f);
        fwrite(&N0,sizeof(int),1,f);
        fwrite(hA,sizeof(float4),N0,f);
    }
    HST void read (FILE* f)
    {   if(hA) free(hA);
        fread (&Z0,sizeof(int),1,f);
        fread (&S0,sizeof(int),1,f);
        fread (&N0,sizeof(int),1,f);
        hA = (float4 *)malloc(N0*sizeof(float4));
        fread(hA,sizeof(float4),N0,f);
        A = hA;
        E0=S0+N0-1;
    }
};

/*
 * Reducing rank of arrays by eliminating last dimension
 *  e.g. A[N][M] = reduce(B[N][M][k=value], where k =0,...K-1 and B=B[N][M][K]
 */

template <class T> HDV kFArray1D<T>& kFArray2D<T>::subMlast(int i1)
{   kFArray1D<T>& M = *(new kFArray1D<T>(S0,E0));
    for(int i0=0;i0<N0;i0++) M.A[i0] = A[(i1-S1)*N0+i0];
    return M;
}
template <class T> HDV kFArray1D<T>& kFArray2D<T>::subMfrst(int i0)
{   kFArray1D<T>& M= *(new kFArray1D<T>(S1,E1));
    for(int i1=0;i1<N1;i1++) M.A[i1] = A[i1*N0+(i0-S0)];
    return M;
}
//////////////////////////////////////////////////////////////////////
template <class T> HDV kFArray1D<T>* kFArray2D<T>::psubMlast(int i1)
{   kFArray1D<T>* P = new kFArray1D<T>(S0,E0);
    for(int i0=0;i0<N0;i0++) P->A[i0] = A[(i1-S1)*N0+i0];
    return P;
}
template <class T> HDV kFArray1D<T>* kFArray2D<T>::psubMfrst(int i0)
{   kFArray1D<T>* P = new kFArray1D<T>(S1,E1);
    for(int i1=0;i1<N1;i1++) P->A[i1] = A[i1*N0+(i0-S0)];
    return P;
}
//////////////////////////////////////////////////////////////////////
template <class T> HDV kCArray1D<T>& kCArray2D<T>::subMlast(int i1)
{   kCArray1D<T>& M = *(new kCArray1D<T>(S0,E0));
    for(int i0=0;i0<N0;i0++) M.A[i0] = A[i0*N1+(i1-S1)];
    return M;
}
template <class T> HDV kCArray1D<T>& kCArray2D<T>::subMfrst(int i0)
{   kCArray1D<T>& M= *(new kCArray1D<T>(S1,E1));
    for(int i1=0;i1<N1;i1++) M.A[i1] = A[(i0-S0)*N1+i1];
    return M;
}

template <class T> HDV kFArray2D<T>& kFArray3D<T>::subMlast(int i2)
{   kFArray2D<T>& M = *(new kFArray2D<T>(S0,E0,S1,E1));
    for(int i1=0;i1<N1;i1++)
    for(int i0=0;i0<N0;i0++) M.A[i1*N0+i0] = A[((i2-S2)*N1+i1)*N0+i0];
    return M;
}
template <class T> HDV kFArray2D<T>& kFArray3D<T>::subMfrst(int i0)
{   kFArray2D<T>& M = *(new kFArray2D<T>(S1,E1,S2,E2));
    for(int i1=0;i1<N1;i1++)
    for(int i2=0;i2<N2;i2++) M.A[i2*N1+i1] = A[(i2*N1+i1)*N0+(i0-S0)];
    return M;
}

template <class T> HDV kCArray2D<T>& kCArray3D<T>::subMlast(int i2)
{   kCArray2D<T>& M = *(new kCArray2D<T>(S0,E0,S1,E1));
    for(int i1=0;i1<N1;i1++)
    for(int i0=0;i0<N0;i0++) M.A[i0*N1+i1] = A[(i0*N1+i1)*N2+(i2-S2)];
    return M;
}
template <class T> HDV kCArray2D<T>& kCArray3D<T>::subMfrst(int i0)
{   kCArray2D<T>& M = *(new kCArray2D<T>(S1,E1,S2,E2));
    for(int i1=0;i1<N1;i1++)
    for(int i2=0;i2<N2;i2++) M.A[i1*N2+i2] = A[((i0-S0)*N1+i1)*N2+i2];
    return M;
}

template <class T> HDV kFArray3D<T>& kFArray4D<T>::subMlast(int i3)
{   kFArray3D<T>& M = *(new kFArray3D<T>(S0,E0,S1,E1,S2,E2));
    for(int i2=0;i2<N2;i2++)
    for(int i1=0;i1<N1;i1++)
    for(int i0=0;i0<N0;i0++) M.A[(i2*N1+i1)*N0+i0] = A[(((i3-S3)*N2+i2)*N1+i1)*N0+i0];
    return M;
}
template <class T> HDV kFArray3D<T>& kFArray4D<T>::subMfrst(int i0)
{   kFArray3D<T>& M = *(new kFArray3D<T>(S1,E1,S2,E2,S3,E3));
    for(int i3=0;i3<N3;i3++)
    for(int i2=0;i2<N2;i2++)
    for(int i1=0;i1<N1;i1++) M.A[(i3*N2+i2)*N1+i1] = A[((i3*N2+i2)*N1+i1)*N0+(i0-S0)];
    return M;
}

template <class T> HDV kCArray3D<T>& kCArray4D<T>::subMlast(int i3)
{   kCArray3D<T>& M = *(new kCArray3D<T>(S0,E0,S1,E1,S2,E2));
    for(int i0=0;i0<N0;i0++)
    for(int i1=0;i1<N1;i1++)
    for(int i2=0;i2<N2;i2++) M.A[(i0*N1+i1)*N2+i2] = A[((i0*N1+i1)*N2+i2)*N3+(i3-S3)];
    return M;
}
template <class T> HDV kCArray3D<T>& kCArray4D<T>::subMfrst(int i0)
{   kCArray3D<T>& M = *(new kCArray3D<T>(S1,E1,S2,E2,S3,E3));
    for(int i1=0;i1<N1;i1++)
    for(int i2=0;i2<N2;i2++)
    for(int i3=0;i3<N3;i3++) M.A[(i1*N2+i2)*N3+i3] = A[(((i0-S0)*N1+i1)*N2+i2)*N3+i3];
    return M;
}

template <class T> HDV kFArray4D<T>& kFArray5D<T>::subMlast(int i4)
{   kFArray4D<T>& M = *(new kFArray4D<T>(S0,E0,S1,E1,S2,E2,S3,E3));
    for(int i3=0;i3<N3;i3++)
    for(int i2=0;i2<N2;i2++)
    for(int i1=0;i1<N1;i1++)
    for(int i0=0;i0<N0;i0++) M.A[((i3*N2+i2)*N1+i1)*N0+i0] = A[((((i4-S4)*N3+i3)*N2+i2)*N1+i1)*N0+i0];
    return M;
}
template <class T> HDV kFArray4D<T>& kFArray5D<T>::subMfrst(int i0)
{   kFArray4D<T>& M = *(new kFArray4D<T>(S1,E1,S2,E2,S3,E3,S4,E4));
    for(int i4=0;i4<N4;i4++)
    for(int i3=0;i3<N3;i3++)
    for(int i2=0;i2<N2;i2++)
    for(int i1=0;i1<N1;i1++) M.A[((i4*N3+i3)*N2+i2)*N1+i1] = A[(((i4*N3+i3)*N2+i2)*N1+i1)*N0+(i0-S0)];
    return M;
}

template <class T> HDV kCArray4D<T>& kCArray5D<T>::subMlast(int i4)
{   kCArray4D<T>& M = *(new kCArray4D<T>(S0,E0,S1,E1,S2,E2,S3,E3));
    for(int i0=0;i0<N0;i0++)
    for(int i1=0;i1<N1;i1++)
    for(int i2=0;i2<N2;i2++)
    for(int i3=0;i3<N3;i3++) M.A[((i0*N1+i1)*N2+i2)*N3+i3] = A[(((i0*N1+i1)*N2+i2)*N3+i3)*N4+(i4-S4)];
    return M;
}
template <class T> HDV kCArray4D<T>& kCArray5D<T>::subMfrst(int i0)
{   kCArray4D<T>& M = *(new kCArray4D<T>(S1,E1,S2,E2,S3,E3,S4,E4));
    for(int i1=0;i1<N1;i1++)
    for(int i2=0;i2<N2;i2++)
    for(int i3=0;i3<N3;i3++)
    for(int i4=0;i4<N4;i4++) M.A[((i1*N2+i2)*N3+i3)*N4+i4] = A[((((i0-S0)*N1+i1)*N2+i2)*N3+i3)*N4+i4];
    return M;
}

template <class T> HDV kFArray5D<T>& kFArray6D<T>::subMlast(int i5)
{   kFArray5D<T>& M = *(new kFArray5D<T>(S0,E0,S1,E1,S2,E2,S3,E3,S4,E4));
    for(int i4=0;i4<N4;i4++)
    for(int i3=0;i3<N3;i3++)
    for(int i2=0;i2<N2;i2++)
    for(int i1=0;i1<N1;i1++)
    for(int i0=0;i0<N0;i0++) M.A[(((i4*N3+i3)*N2+i2)*N1+i1)*N0+i0] = A[(((((i5-S5)*N4+i4)*N3+i3)*N2+i2)*N1+i1)*N0+i0];
    return M;
}
template <class T> HDV kFArray5D<T>& kFArray6D<T>::subMfrst(int i0)
{   kFArray5D<T>& M = *(new kFArray5D<T>(S1,E1,S2,E2,S3,E3,S4,E4,S5,E5));
    for(int i5=0;i5<N5;i5++)
    for(int i4=0;i4<N4;i4++)
    for(int i3=0;i3<N3;i3++)
    for(int i2=0;i2<N2;i2++)
    for(int i1=0;i1<N1;i1++) M.A[(((i5*N4+i4)*N3+i3)*N2+i2)*N1+i1] = A[((((i5*N4+i4)*N3+i3)*N2+i2)*N1+i1)*N0+(i0-S0)];
    return M;
}

template <class T> HDV kCArray5D<T>& kCArray6D<T>::subMlast(int i5)
{   kCArray5D<T>& M = *(new kCArray5D<T>(S0,E0,S1,E1,S2,E2,S3,E3,S4,E4));
    for(int i0=0;i0<N0;i0++)
    for(int i1=0;i1<N1;i1++)
    for(int i2=0;i2<N2;i2++)
    for(int i3=0;i3<N3;i3++)
    for(int i4=0;i4<N4;i4++) M.A[(((i0*N1+i1)*N2+i2)*N3+i3)*N4+i4] = A[((((i0*N1+i1)*N2+i2)*N3+i3)*N4+i4)*N5+(i5-S5)];
    return M;
}
template <class T> HDV kCArray5D<T>& kCArray6D<T>::subMfrst(int i0)
{   kCArray5D<T>& M = *(new kCArray5D<T>(S1,E1,S2,E2,S3,E3,S4,E4,S5,E5));
    for(int i1=0;i1<N1;i1++)
    for(int i2=0;i2<N2;i2++)
    for(int i3=0;i3<N3;i3++)
    for(int i4=0;i4<N4;i4++)
    for(int i5=0;i5<N5;i5++) M.A[(((i1*N2+i2)*N3+i3)*N4+i4)*N5+i5] = A[(((((i0-S0)*N1+i1)*N2+i2)*N3+i3)*N4+i4)*N5+i5];
    return M;
}

template <class T> HDV kFArray6D<T>& kFArray7D<T>::subMlast(int i6)
{   kFArray6D<T>& M = *(new kFArray6D<T>(S0,E0,S1,E1,S2,E2,S3,E3,S4,E4,S5,E5));
    for(int i5=0;i5<N5;i5++)
    for(int i4=0;i4<N4;i4++)
    for(int i3=0;i3<N3;i3++)
    for(int i2=0;i2<N2;i2++)
    for(int i1=0;i1<N1;i1++)
    for(int i0=0;i0<N0;i0++) M.A[((((i5*N4+i4)*N3+i3)*N2+i2)*N1+i1)*N0+i0] = A[((((((i6-S6)*N5+i5)*N4+i4)*N3+i3)*N2+i2)*N1+i1)*N0+i0];
    return M;
}
template <class T> HDV kFArray6D<T>& kFArray7D<T>::subMfrst(int i0)
{   kFArray6D<T>& M = *(new kFArray6D<T>(S1,E1,S2,E2,S3,E3,S4,E4,S5,E5,S6,E6));
    for(int i6=0;i6<N6;i6++)
    for(int i5=0;i5<N5;i5++)
    for(int i4=0;i4<N4;i4++)
    for(int i3=0;i3<N3;i3++)
    for(int i2=0;i2<N2;i2++)
    for(int i1=0;i1<N1;i1++) M.A[((((i6*N5+i5)*N4+i4)*N3+i3)*N2+i2)*N1+i1] = A[(((((i6*N5+i5)*N4+i4)*N3+i3)*N2+i2)*N1+i1)*N0+(i0-S0)];
    return M;
}

template <class T> HDV kCArray6D<T>& kCArray7D<T>::subMlast(int i6)
{   kCArray6D<T>& M = *(new kCArray6D<T>(S0,E0,S1,E1,S2,E2,S3,E3,S4,E4,S5,E5));
    for(int i0=0;i0<N0;i0++)
    for(int i1=0;i1<N1;i1++)
    for(int i2=0;i2<N2;i2++)
    for(int i3=0;i3<N3;i3++)
    for(int i4=0;i4<N4;i4++)
    for(int i5=0;i5<N5;i5++) M.A[((((i0*N1+i1)*N2+i2)*N3+i3)*N4+i4)*N5+i5] = A[(((((i0*N1+i1)*N2+i2)*N3+i3)*N4+i4)*N5+i5)*N6+(i6-S6)];
    return M;
}
template <class T> HDV kCArray6D<T>& kCArray7D<T>::subMfrst(int i0)
{   kCArray6D<T>& M = *(new kCArray6D<T>(S1,E1,S2,E2,S3,E3,S4,E4,S5,E5,S6,E6));
    for(int i1=0;i1<N1;i1++)
    for(int i2=0;i2<N2;i2++)
    for(int i3=0;i3<N3;i3++)
    for(int i4=0;i4<N4;i4++)
    for(int i5=0;i5<N5;i5++)
    for(int i6=0;i6<N6;i6++) M.A[((((i1*N2+i2)*N3+i3)*N4+i4)*N5+i5)*N6+i6] = A[((((((i0-S0)*N1+i1)*N2+i2)*N3+i3)*N4+i4)*N5+i5)*N6+i6];
    return M;
}

/*
 * Matrix Transfer Operations - to/from device
 */
#if NVCC
//1D matrices
template <class type> HST kMatrix1D<type>* copyMatrixHostToDevice(kMatrix1D<type>* hM)
{
    int Nelmts = hM->size();
    kMatrix1D<type>* dM;
    type* hA = hM->hA;
    type* dA;
    cudaMalloc((void **) &dA,Nelmts*sizeof(type));
    cudaMemcpy(dA, hA, Nelmts*sizeof(type), cudaMemcpyHostToDevice);
    hM->dA= dA;
    hM->hA= hA;
    hM->A = dA;
    cudaMalloc((void **) &dM,sizeof(*hM));
    cudaMemcpy(dM, hM, sizeof(*hM), cudaMemcpyHostToDevice);
    hM->A = hA;
    return dM;
}
template <class type> HST kMatrix1D<type>* copyMatrixDeviceToHost(kMatrix1D<type>* dM)
{
    kMatrix1D<type>* hM= new kMatrix1D<type>();
    cudaMemcpy(hM, dM, sizeof(*hM), cudaMemcpyDeviceToHost);
    int Nelmts = hM->size();
    type* dA = hM->A;
    type* hA = (type *)malloc(Nelmts*sizeof(type));
    cudaMemcpy(hA, dA, Nelmts*sizeof(type), cudaMemcpyDeviceToHost);
    hM->hA= hA;
    hM->dA= dA;
    hM->A = hA;
    return hM;
}
//2D matrices
template <class type> HST kMatrix2D<type>* copyMatrixHostToDevice(kMatrix2D<type>* hM)
{
    int Nelmts = hM->size();
    kMatrix2D<type>* dM;
    type* hA = hM->hA;
    type* dA;
    cudaMalloc((void **) &dA,Nelmts*sizeof(type));
    cudaMemcpy(dA, hA, Nelmts*sizeof(type), cudaMemcpyHostToDevice);
    hM->dA= dA;
    hM->hA= hA;
    hM->A = dA;
    cudaMalloc((void **) &dM,sizeof(*hM));
    cudaMemcpy(dM, hM, sizeof(*hM), cudaMemcpyHostToDevice);
    hM->A = hA;
    return dM;
}
template <class type> HST kMatrix2D<type>* copyMatrixDeviceToHost(kMatrix2D<type>* dM)
{
    kMatrix2D<type>* hM= new kMatrix2D<type>();
    cudaMemcpy(hM, dM, sizeof(*hM), cudaMemcpyDeviceToHost);;
    int Nelmts = hM->size();
    type* dA = hM->A;
    type* hA = (type *)malloc(Nelmts*sizeof(type));
    cudaMemcpy(hA, dA, Nelmts*sizeof(type), cudaMemcpyDeviceToHost);
    hM->hA= hA;
    hM->dA= dA;
    hM->A = hA;
    return hM;
}
//3D matrices
template <class type> HST kMatrix3D<type>* copyMatrixHostToDevice(kMatrix3D<type>* hM)
{
    int Nelmts = hM->size();
    kMatrix3D<type>* dM;
    type* hA = hM->hA;
    type* dA;
    cudaMalloc((void **) &dA,Nelmts*sizeof(type));
    cudaMemcpy(dA, hA, Nelmts*sizeof(type), cudaMemcpyHostToDevice);
    hM->dA= dA;
    hM->hA= hA;
    hM->A = dA;
    cudaMalloc((void **) &dM,sizeof(*hM));
    cudaMemcpy(dM, hM, sizeof(*hM), cudaMemcpyHostToDevice);
    hM->A = hA;
    return dM;
}
template <class type> HST kMatrix3D<type>* copyMatrixDeviceToHost(kMatrix3D<type>* dM)
{
    kMatrix3D<type>* hM= new kMatrix3D<type>();
    cudaMemcpy(hM, dM, sizeof(*hM), cudaMemcpyDeviceToHost);
    int Nelmts = hM->size();
    type* dA = hM->A;
    type* hA = (type *)malloc(Nelmts*sizeof(type));
    cudaMemcpy(hA, dA, Nelmts*sizeof(type), cudaMemcpyDeviceToHost);
    hM->hA= hA;
    hM->dA= dA;
    hM->A = hA;
    return hM;
}
//4D matrices
template <class type> HST kMatrix4D<type>* copyMatrixHostToDevice(kMatrix4D<type>* hM)
{
    int Nelmts = hM->size();
    kMatrix4D<type>* dM;
    type* hA = hM->hA;
    type* dA;
    cudaMalloc((void **) &dA,Nelmts*sizeof(type));
    cudaMemcpy(dA, hA, Nelmts*sizeof(type), cudaMemcpyHostToDevice);
    hM->dA= dA;
    hM->hA= hA;
    hM->A = dA;
    cudaMalloc((void **) &dM,sizeof(*hM));
    cudaMemcpy(dM, hM, sizeof(*hM), cudaMemcpyHostToDevice);
    hM->A = hA;
    return dM;
}
template <class type> HST kMatrix4D<type>* copyMatrixDeviceToHost(kMatrix4D<type>* dM)
{
    kMatrix4D<type>* hM= new kMatrix4D<type>();
    cudaMemcpy(hM, dM, sizeof(*hM), cudaMemcpyDeviceToHost);
    int Nelmts = hM->size();
    type* dA = hM->A;
    type* hA = (type *)malloc(Nelmts*sizeof(type));
    cudaMemcpy(hA, dA, Nelmts*sizeof(type), cudaMemcpyDeviceToHost);
    hM->hA= hA;
    hM->dA= dA;
    hM->A = hA;
    return hM;
}
//5D matrices
template <class type> HST kMatrix5D<type>* copyMatrixHostToDevice(kMatrix5D<type>* hM)
{
    int Nelmts = hM->size();
    kMatrix5D<type>* dM;
    type* hA = hM->hA;
    type* dA;
    cudaMalloc((void **) &dA,Nelmts*sizeof(type));
    cudaMemcpy(dA, hA, Nelmts*sizeof(type), cudaMemcpyHostToDevice);
    hM->dA= dA;
    hM->hA= hA;
    hM->A = dA;
    cudaMalloc((void **) &dM,sizeof(*hM));
    cudaMemcpy(dM, hM, sizeof(*hM), cudaMemcpyHostToDevice);
    hM->A = hA;
    return dM;
}
template <class type> HST kMatrix5D<type>* copyMatrixDeviceToHost(kMatrix5D<type>* dM)
{
    kMatrix5D<type>* hM= new kMatrix5D<type>();
    cudaMemcpy(hM, dM, sizeof(*hM), cudaMemcpyDeviceToHost);
    int Nelmts = hM->size();
    type* dA = hM->A;
    type* hA = (type *)malloc(Nelmts*sizeof(type));
    cudaMemcpy(hA, dA, Nelmts*sizeof(type), cudaMemcpyDeviceToHost);
    hM->hA = hA;
    hM->dA = dA;
    hM->A  = hA;
    return hM;
}
//6D matrices
template <class type> HST kMatrix6D<type>* copyMatrixHostToDevice(kMatrix6D<type>* hM)
{
    int Nelmts = hM->size();
    kMatrix6D<type>* dM;
    type* hA = hM->hA;
    type* dA;
    cudaMalloc((void **) &dA,Nelmts*sizeof(type));
    cudaMemcpy(dA, hA, Nelmts*sizeof(type), cudaMemcpyHostToDevice);
    hM->dA= dA;
    hM->hA= hA;
    hM->A = dA;
    cudaMalloc((void **) &dM,sizeof(*hM));
    cudaMemcpy(dM, hM, sizeof(*hM), cudaMemcpyHostToDevice);
    hM->A = hA;
    return dM;
}
template <class type> HST kMatrix6D<type>* copyMatrixDeviceToHost(kMatrix6D<type>* dM)
{
    kMatrix6D<type>* hM= new kMatrix6D<type>();
    cudaMemcpy(hM, dM, sizeof(*hM), cudaMemcpyDeviceToHost);
    int Nelmts = hM->size();
    type* dA = hM->A;
    type* hA = (type *)malloc(Nelmts*sizeof(type));
    cudaMemcpy(hA, dA, Nelmts*sizeof(type), cudaMemcpyDeviceToHost);
    hM->hA= hA;
    hM->dA= dA;
    hM->A = hA;
    return hM;
}
//7D matrices
template <class type> HST kMatrix7D<type>* copyMatrixHostToDevice(kMatrix7D<type>* hM)
{
    int Nelmts = hM->size();
    kMatrix7D<type>* dM;
    type* hA = hM->hA;
    type* dA;
    cudaMalloc((void **) &dA,Nelmts*sizeof(type));
    cudaMemcpy(dA, hA, Nelmts*sizeof(type), cudaMemcpyHostToDevice);
    hM->dA= dA;
    hM->hA= hA;
    hM->A = dA;
    cudaMalloc((void **) &dM,sizeof(*hM));
    cudaMemcpy(dM, hM, sizeof(*hM), cudaMemcpyHostToDevice);
    hM->A = hA;
    return dM;
}
template <class type> HST kMatrix7D<type>* copyMatrixDeviceToHost(kMatrix7D<type>* dM)
{
    kMatrix7D<type>* hM= new kMatrix7D<type>();
    cudaMemcpy(hM, dM, sizeof(*hM), cudaMemcpyDeviceToHost);
    int Nelmts = hM->size();
    type* dA = hM->A;
    type* hA = (type *)malloc(Nelmts*sizeof(type));
    cudaMemcpy(hA, dA, Nelmts*sizeof(type), cudaMemcpyDeviceToHost);
    hM->hA= hA;
    hM->dA= dA;
    hM->A = hA;
    return hM;
}
#endif

/*
*  library functions: dot products (contraction of the last indices)
*    ex: A[i][j] = dot(B[i][j][k],C[k])
*    if same dimensions -> dot product: a = dot(A[i][j],B[i][j])
*/

template<class T> HDV T dot(const kFArray1D<T>& L, const kFArray1D<T>& R)
{   T prod = (T)(0.);
    for(int i=0;i<R.N0;i++)  prod += (L.A[i]*R.A[i]);
    return prod;
}
template<class T> HDV T dot(const kCArray1D<T>& L, const kCArray1D<T>& R)
{   T prod = (T)(0.);
    for(int i=0;i<R.N0;i++)  prod += (L.A[i]*R.A[i]);
    return prod;
}

template<class T> HDV T dot(const kFArray2D<T>& L, const kFArray2D<T>& R)
{   T prod = (T)(0.);
    for(int i=0;i<R.N0*R.N1;i++) prod += L.A[i]*R.A[i];
    return prod;
}
template<class T> HDV T dot(const kCArray2D<T>& L, const kCArray2D<T>& R)
{   T prod = (T)(0.);
    for(int i=0;i<R.N0*R.N1;i++) prod += L.A[i]*R.A[i];
    return prod;
}

template<class T> HDV T dot(const kFArray3D<T>& L, const kFArray3D<T>& R)
{   T prod = (T)(0.);
    for(int i=0;i<R.N0*R.N1*R.N2;i++) prod += L.A[i]*R.A[i];
    return prod;
}
template<class T> HDV T dot(const kCArray3D<T>& L, const kCArray3D<T>& R)
{   T prod = (T)(0.);
    for(int i=0;i<R.N0*R.N1*R.N2;i++) prod += L.A[i]*R.A[i];
    return prod;
}

template<class T> HDV T dot(const kFArray4D<T>& L, const kFArray4D<T>& R)
{   T prod = (T)(0.);
    for(int i=0;i<R.N0*R.N1*R.N2*R.N3;i++) prod += L.A[i]*R.A[i];
    return prod;
}
template<class T> HDV T dot(const kCArray4D<T>& L, const kCArray4D<T>& R)
{   T prod = (T)(0.);
    for(int i=0;i<R.N0*R.N1*R.N2*R.N3;i++) prod += L.A[i]*R.A[i];
    return prod;
}

template<class T> HDV T dot(const kFArray5D<T>& L, const kFArray5D<T>& R)
{   T prod = (T)(0.);
    for(int i=0;i<R.N0*R.N1*R.N2*R.N3*R.N4;i++) prod += L.A[i]*R.A[i];
    return prod;
}
template<class T> HDV T dot(const kCArray5D<T>& L, const kCArray5D<T>& R)
{   T prod = (T)(0.);
    for(int i=0;i<R.N0*R.N1*R.N2*R.N3*R.N4;i++) prod += L.A[i]*R.A[i];
    return prod;
}

template<class T> HDV T dot(const kFArray6D<T>& L, const kFArray6D<T>& R)
{   T prod = (T)(0.);
    for(int i=0;i<R.N0*R.N1*R.N2*R.N3*R.N4*R.N5;i++) prod += L.A[i]*R.A[i];
    return prod;
}
template<class T> HDV T dot(const kCArray6D<T>& L, const kCArray6D<T>& R)
{   T prod = (T)(0.);
    for(int i=0;i<R.N0*R.N1*R.N2*R.N3*R.N4*R.N5;i++) prod += L.A[i]*R.A[i];
    return prod;
}

template<class T> HDV T dot(const kFArray7D<T>& L, const kFArray7D<T>& R)
{   T prod = (T)(0.);
    for(int i=0;i<R.N0*R.N1*R.N2*R.N3*R.N4*R.N5*R.N6;i++) prod += L.A[i]*R.A[i];
    return prod;
}
template<class T> HDV T dot(const kCArray7D<T>& L, const kCArray7D<T>& R)
{   T prod = (T)(0.);
    for(int i=0;i<R.N0*R.N1*R.N2*R.N3*R.N4*R.N5*R.N6;i++) prod += L.A[i]*R.A[i];
    return prod;
}

// dot products returning 1D arrays
template<class T> HDV kFArray1D<T>* dot(const kFArray2D<T>& L, const kFArray1D<T>& R)
{   kFArray1D<T>* P = new kFArray1D<T>(L.S0,L.E0);
    for(int i0=0;i0<L.N0;i0++)
    {   P->A[i0]=(T) 0.;
        for(int i1=0;i1<L.N1;i1++) P->A[i0] += (L.A[i1*L.N0+i0]*R.A[i1]);
    }
    return P;
}
template<class T> HDV kCArray1D<T>* dot(const kCArray2D<T>& L, const kCArray1D<T>& R)
{   kCArray1D<T>* P = new kCArray1D<T>(L.S0,L.E0);
    for(int i0=0;i0<L.N0;i0++)
    {   P->A[i0]=(T) 0.;
        for(int i1=0;i1<L.N1;i1++) P->A[i0] += (L.A[i0*L.N1+i1]*R.A[i1]);
    }
    return P;
}

template<class T> HDV kFArray1D<T>* dot(const kFArray3D<T>& L, const kFArray2D<T>& R)
{   kFArray1D<T>* P = new kFArray1D<T>(L.S0,L.E0);
    for(int i0=0;i0<L.N0;i0++)
    {   P->A[i0]=(T) 0.;  // assumes R.N0 = L.N1 and R.N1 = L.N2
        for(int i1=0;i1<L.N1;i1++)
        for(int i2=0;i2<L.N2;i2++) P->A[i0] += (L.A[(i2*L.N1+i1)*L.N0+i0]*R.A[i2*L.N1+i1]);
    }
    return P;
}
template<class T> HDV kCArray1D<T>* dot(const kCArray3D<T>& L, const kCArray2D<T>& R)
{   kCArray1D<T>* P = new kCArray1D<T>(L.S0,L.E0);
    for(int i0=0;i0<L.N0;i0++)
    {   P->A[i0]=(T) 0.;  // assumes R.N0 = L.N1 and R.N1 = L.N2
        for(int i1=0;i1<L.N1;i1++)
        for(int i2=0;i2<L.N2;i2++) P->A[i0] += (L.A[(i0*L.N1+i1)*L.N2+i2]*R.A[i1*L.N2+i2]);
    }
    return P;
}

template<class T> HDV kFArray1D<T>* dot(const kFArray4D<T>& L, const kFArray3D<T>& R)
{   kFArray1D<T>* P = new kFArray1D<T>(L.S0,L.E0);
    for(int i0=0;i0<L.N0;i0++)
    {   P->A[i0]=(T) 0.;  // assumes R.N0 = L.N1, R.N1 = L.N2 and R.N2 = L.N3
        for(int i1=0;i1<L.N1;i1++)
        for(int i2=0;i2<L.N2;i2++)
        for(int i3=0;i3<L.N3;i3++) P->A[i0] += (L.A[((i3*L.N2+i2)*L.N1+i1)*L.N0+i0]*R.A[(i3*L.N2+i2)*L.N1+i1]);
    }
    return P;
}
template<class T> HDV kCArray1D<T>* dot(const kCArray4D<T>& L, const kCArray3D<T>& R)
{   kCArray1D<T>* P = new kCArray1D<T>(L.S0,L.E0);
    for(int i0=0;i0<L.N0;i0++)
    {   P->A[i0]=(T) 0.;  // assumes R.N0 = L.N1, R.N1 = L.N2 and R.N2 = L.N3
        for(int i1=0;i1<L.N1;i1++)
        for(int i2=0;i2<L.N2;i2++)
        for(int i3=0;i3<L.N3;i3++) P->A[i0] += (L.A[((i0*L.N1+i1)*L.N2+i2)*L.N3+i3]*R.A[(i1*L.N2+i2)*L.N3+i3]);
    }
    return P;
}

template<class T> HDV kFArray1D<T>* dot(const kFArray5D<T>& L, const kFArray4D<T>& R)
{   kFArray1D<T>* P = new kFArray1D<T>(L.S0,L.E0);
    for(int i0=0;i0<L.N0;i0++)
    {   P->A[i0]=(T) 0.;  // assumes R.N0 = L.N1, R.N1 = L.N2, R.N2 = L.N3 and R.N3 = L.N4
        for(int i1=0;i1<L.N1;i1++)
        for(int i2=0;i2<L.N2;i2++)
        for(int i3=0;i3<L.N3;i3++)
        for(int i4=0;i4<L.N4;i4++) P->A[i0] += (L.A[(((i4*L.N3+i3)*L.N2+i2)*L.N1+i1)*L.N0+i0]*R.A[((i4*L.N3+i3)*L.N2+i2)*L.N1+i1]);
    }
    return P;
}
template<class T> HDV kCArray1D<T>* dot(const kCArray5D<T>& L, const kCArray4D<T>& R)
{   kCArray1D<T>* P = new kCArray1D<T>(L.S0,L.E0);
    for(int i0=0;i0<L.N0;i0++)
    {   P->A[i0]=(T) 0.;  // assumes R.N0 = L.N1, R.N1 = L.N2, R.N2 = L.N3 and R.N3 = L.N4
        for(int i1=0;i1<L.N1;i1++)
        for(int i2=0;i2<L.N2;i2++)
        for(int i3=0;i3<L.N3;i3++)
        for(int i4=0;i4<L.N4;i4++) P->A[i0] += (L.A[(((i0*L.N1+i1)*L.N2+i2)*L.N3+i3)*L.N4+i4]*R.A[((i1*L.N2+i2)*L.N3+i3)*L.N4+i4]);
    }
    return P;
}

template<class T> HDV kFArray1D<T>* dot(const kFArray6D<T>& L, const kFArray5D<T>& R)
{   kFArray1D<T>* P = new kFArray1D<T>(L.S0,L.E0);
    for(int i0=0;i0<L.N0;i0++)
    {   P->A[i0]=(T) 0.;  // assumes R.N0 = L.N1, R.N1 = L.N2, R.N2 = L.N3, R.N3 = L.N4 and R.N4 = L.N5
        for(int i1=0;i1<L.N1;i1++)
        for(int i2=0;i2<L.N2;i2++)
        for(int i3=0;i3<L.N3;i3++)
        for(int i4=0;i4<L.N4;i4++)
        for(int i5=0;i5<L.N5;i5++) P->A[i0] += (L.A[((((i5*L.N4+i4)*L.N3+i3)*L.N2+i2)*L.N1+i1)*L.N0+i0]*R.A[(((i5*L.N4+i4)*L.N3+i3)*L.N2+i2)*L.N1+i1]);
    }
    return P;
}
template<class T> HDV kCArray1D<T>* dot(const kCArray6D<T>& L, const kCArray5D<T>& R)
{   kCArray1D<T>* P = new kCArray1D<T>(L.S0,L.E0);
    for(int i0=0;i0<L.N0;i0++)
    {   P->A[i0]=(T) 0.;  // assumes R.N0 = L.N1, R.N1 = L.N2, R.N2 = L.N3, R.N3 = L.N4 and R.N4 = L.N5
        for(int i1=0;i1<L.N1;i1++)
        for(int i2=0;i2<L.N2;i2++)
        for(int i3=0;i3<L.N3;i3++)
        for(int i4=0;i4<L.N4;i4++)
        for(int i5=0;i5<L.N5;i5++) P->A[i0] += (L.A[((((i0*L.N1+i1)*L.N2+i2)*L.N3+i3)*L.N4+i4)*L.N5+i5]*R.A[(((i1*L.N2+i2)*L.N3+i3)*L.N4+i4)*L.N5+i5]);
    }
    return P;
}

template<class T> HDV kFArray1D<T>* dot(const kFArray7D<T>& L, const kFArray6D<T>& R)
{   kFArray1D<T>* P = new kFArray1D<T>(L.S0,L.E0);
    for(int i0=0;i0<L.N0;i0++)
    {   P->A[i0]=(T) 0.;  // assumes R.N0 = L.N1, R.N1 = L.N2, R.N2 = L.N3, R.N3 = L.N4, R.N4 = L.N5 and R.N5=L.N6
        for(int i1=0;i1<L.N1;i1++)
        for(int i2=0;i2<L.N2;i2++)
        for(int i3=0;i3<L.N3;i3++)
        for(int i4=0;i4<L.N4;i4++)
        for(int i5=0;i5<L.N5;i5++)
        for(int i6=0;i6<L.N6;i6++) P->A[i0] += (L.A[(((((i6*L.N5+i5)*L.N4+i4)*L.N3+i3)*L.N2+i2)*L.N1+i1)*L.N0+i0]*R.A[((((i6*L.N5+i5)*L.N4+i4)*L.N3+i3)*L.N2+i2)*L.N1+i1]);
    }
    return P;
}
template<class T> HDV kCArray1D<T>* dot(const kCArray7D<T>& L, const kCArray6D<T>& R)
{   kFArray1D<T>* P = new kFArray1D<T>(L.S0,L.E0);
    for(int i0=0;i0<L.N0;i0++)
    {   P->A[i0]=(T) 0.;  // assumes R.N0 = L.N1, R.N1 = L.N2, R.N2 = L.N3, R.N3 = L.N4, R.N4 = L.N5 and R.N5=L.N6
        for(int i1=0;i1<L.N1;i1++)
        for(int i2=0;i2<L.N2;i2++)
        for(int i3=0;i3<L.N3;i3++)
        for(int i4=0;i4<L.N4;i4++)
        for(int i5=0;i5<L.N5;i5++)
        for(int i6=0;i6<L.N6;i6++) P->A[i0] += (L.A[(((((i0*L.N1+i1)*L.N2+i2)*L.N3+i3)*L.N4+i4)*L.N5+i5)*L.N6+i6]*R.A[((((i1*L.N2+i2)*L.N3+i3)*L.N4+i4)*L.N5+i5)*L.N6+i6]);
    }
    return P;
}

// dot products returning 2D matrix
template<class T> HDV kFArray2D<T>* dot(const kFArray3D<T>& L, const kFArray1D<T>& R)
{   kFArray2D<T>* P = new kFArray2D<T>(L.S0,L.E0,L.S1,L.E1);
    for(int i0=0;i0<L.N0;i0++)
    for(int i1=0;i1<L.N1;i1++)
    {   P->A[i1*L.N0+i0]=(T) 0.; // assumes R.N0=L.N2
        for(int i2=0;i2<R.N0;i2++) P->A[i1*L.N0+i0] += (L.A[(i2*L.N1+i1)*L.N0+i0]*R.A[i2]);
    }
    return P;
}
template<class T> HDV kCArray2D<T>* dot(const kCArray3D<T>& L, const kCArray1D<T>& R)
{   kCArray2D<T>* P = new kCArray2D<T>(L.S0,L.E0,L.S1,L.E1);
    for(int i0=0;i0<L.N0;i0++)
    for(int i1=0;i1<L.N1;i1++)
    {   P->A[i0*L.N1+i1]=(T) 0.; // assumes R.N0=L.N2
        for(int i2=0;i2<L.N2;i2++) P->A[i0*L.N1+i1] += (L.A[(i0*L.N1+i1)*L.N2+i2]*R.A[i2]);
    }
    return P;
}

template<class T> HDV kFArray2D<T>* dot(const kFArray4D<T>& L, const kFArray2D<T>& R)
{   kFArray2D<T>* P = new kFArray2D<T>(L.S0,L.E0,L.S1,L.E1);
    for(int i0=0;i0<L.N0;i0++)
    for(int i1=0;i1<L.N1;i1++)
    {   P->A[i1*L.N0+i0]=(T) 0.; // assumes R.N0=L.N2, R.N1=L.N3
        for(int i2=0;i2<L.N2;i2++)
        for(int i3=0;i3<L.N3;i3++) P->A[i1*L.N0+i0] += (L.A[((i3*L.N2+i2)*L.N1+i1)*L.N0+i0]*R.A[i3*L.N2+i2]);
    }
    return P;
}
template<class T> HDV kCArray2D<T>* dot(const kCArray4D<T>& L, const kCArray2D<T>& R)
{   kCArray2D<T>* P = new kCArray2D<T>(L.S0,L.E0,L.S1,L.E1);
    for(int i0=0;i0<L.N0;i0++)
    for(int i1=0;i1<L.N1;i1++)
    {   P->A[i0*L.N1+i1]=(T) 0.; // assumes R.N0=L.N2, R.N1=L.N3
        for(int i2=0;i2<L.N2;i2++)
        for(int i3=0;i3<L.N3;i3++) P->A[i0*L.N1+i1] += (L.A[((i0*L.N1+i1)*L.N2+i2)*L.N3+i3]*R.A[i2*L.N3+i3]);
    }
    return P;
}

template<class T> HDV kFArray2D<T>* dot(const kFArray5D<T>& L, const kFArray3D<T>& R)
{   kFArray2D<T>* P = new kFArray2D<T>(L.S0,L.E0,L.S1,L.E1);
    for(int i0=0;i0<L.N0;i0++)
    for(int i1=0;i1<L.N1;i1++)
    {   P->A[i1*L.N0+i0]=(T) 0.; // assumes R.N0=L.N2, R.N1=L.N3, R.N2=L.N4
        for(int i2=0;i2<L.N2;i2++)
        for(int i3=0;i3<L.N3;i3++)
        for(int i4=0;i4<L.N4;i4++) P->A[i1*L.N0+i0] += (L.A[(((i4*L.N3+i3)*L.N2+i2)*L.N1+i1)*L.N0+i0]*R.A[(i4*L.N3+i3)*L.N2+i2]);
    }
    return P;
}
template<class T> HDV kCArray2D<T>* dot(const kCArray5D<T>& L, const kCArray3D<T>& R)
{   kCArray2D<T>* P = new kCArray2D<T>(L.S0,L.E0,L.S1,L.E1);
    for(int i0=0;i0<L.N0;i0++)
    for(int i1=0;i1<L.N1;i1++)
    {   P->A[i0*L.N1+i1]=(T) 0.; // assumes R.N0=L.N2, R.N1=L.N3, R.N2=L.N4
        for(int i2=0;i2<L.N2;i2++)
        for(int i3=0;i3<L.N3;i3++)
        for(int i4=0;i4<L.N4;i4++) P->A[i0*L.N1+i1] += (L.A[(((i0*L.N1+i1)*L.N2+i2)*L.N3+i3)*L.N4+i4]*R.A[(i2*L.N3+i3)*L.N4+i4]);
    }
    return P;
}

template<class T> HDV kFArray2D<T>* dot(const kFArray6D<T>& L, const kFArray4D<T>& R)
{   kFArray2D<T>* P = new kFArray2D<T>(L.S0,L.E0,L.S1,L.E1);
    for(int i0=0;i0<L.N0;i0++)
    for(int i1=0;i1<L.N1;i1++)
    {   P->A[i1*L.N0+i0]=(T) 0.; // assumes R.N0=L.N2, R.N1=L.N3, R.N2=L.N4, R.N3=L.N5
        for(int i2=0;i2<L.N2;i2++)
        for(int i3=0;i3<L.N3;i3++)
        for(int i4=0;i4<L.N4;i4++)
        for(int i5=0;i5<L.N5;i5++) P->A[i1*L.N0+i0] += (L.A[((((i5*L.N4+i4)*L.N3+i3)*L.N2+i2)*L.N1+i1)*L.N0+i0]*R.A[((i5*L.N4+i4)*L.N3+i3)*L.N2+i2]);
    }
    return P;
}
template<class T> HDV kCArray2D<T>* dot(const kCArray6D<T>& L, const kCArray4D<T>& R)
{   kCArray2D<T>* P = new kCArray2D<T>(L.S0,L.E0,L.S1,L.E1);
    for(int i0=0;i0<L.N0;i0++)
    for(int i1=0;i1<L.N1;i1++)
    {   P->A[i0*L.N1+i1]=(T) 0.; // assumes R.N0=L.N2, R.N1=L.N3, R.N2=L.N4, R.N3=L.N5
        for(int i2=0;i2<L.N2;i2++)
        for(int i3=0;i3<L.N3;i3++)
        for(int i4=0;i4<L.N4;i4++)
        for(int i5=0;i5<L.N5;i5++) P->A[i0*L.N1+i1] += (L.A[((((i0*L.N1+i1)*L.N2+i2)*L.N3+i3)*L.N4+i4)*L.N5+i5]*R.A[((i2*L.N3+i3)*L.N4+i4)*L.N5+i5]);
    }
    return P;
}

template<class T> HDV kFArray2D<T>* dot(const kFArray7D<T>& L, const kFArray5D<T>& R)
{   kFArray2D<T>* P = new kFArray2D<T>(L.S0,L.E0,L.S1,L.E1);
    for(int i0=0;i0<L.N0;i0++)
    for(int i1=0;i1<L.N1;i1++)
    {   P->A[i1*L.N0+i0]=(T) 0.; // assumes R.N0=L.N2, R.N1=L.N3, R.N2=L.N4, R.N3=L.N5, R.N4=L.N6
        for(int i2=0;i2<L.N6;i2++)
        for(int i3=0;i3<L.N5;i3++)
        for(int i4=0;i4<L.N4;i4++)
        for(int i5=0;i5<L.N3;i5++)
        for(int i6=0;i6<L.N2;i6++) P->A[i1*L.N0+i0] += (L.A[(((((i6*L.N5+i5)*L.N4+i4)*L.N3+i3)*L.N2+i2)*L.N1+i1)*L.N0+i0]*R.A[(((i6*L.N5+i5)*L.N4+i4)*L.N3+i3)*L.N2+i2]);
    }
    return P;
}
template<class T> HDV kCArray2D<T>* dot(const kCArray7D<T>& L, const kCArray5D<T>& R)
{   kCArray2D<T>* P = new kCArray2D<T>(L.S0,L.E0,L.S1,L.E1);
    for(int i1=0;i1<L.N1;i1++)
    for(int i0=0;i0<L.N0;i0++)
    {   P->A[i0*L.N1+i1]=(T) 0.; // assumes R.N0=L.N2, R.N1=L.N3, R.N2=L.N4, R.N3=L.N5, R.N4=L.N6
        for(int i2=0;i2<L.N6;i2++)
        for(int i3=0;i3<L.N5;i3++)
        for(int i4=0;i4<L.N4;i4++)
        for(int i5=0;i5<L.N3;i5++)
        for(int i6=0;i6<L.N2;i6++) P->A[i0*L.N0+i1] += (L.A[(((((i0*L.N1+i1)*L.N2+i2)*L.N3+i3)*L.N4+i4)*L.N5+i5)*L.N6+i6]*R.A[(((i2*L.N3+i3)*L.N4+i4)*L.N5+i5)*L.N6+i6]);
    }
    return P;
}

// dot products returning 3D matrix
template<class T> HDV kFArray3D<T>* dot(const kFArray4D<T>& L, const kFArray1D<T>& R)
{   kFArray3D<T>* P = new kFArray3D<T>(L.S0,L.E0,L.S1,L.E1,L.S2,L.E2);
    for(int i0=0;i0<L.N0;i0++)
    for(int i1=0;i1<L.N1;i1++)
    for(int i2=0;i2<L.N2;i2++)
    {   P->A[(i2*L.N1+i1)*L.N0+i0]=(T) 0.; // assumes R.N0=L.N3
        for(int i3=0;i3<L.N3;i3++) P->A[(i2*L.N1+i1)*L.N0+i0] += (L.A[((i3*L.N2+i2)*L.N1+i1)*L.N0+i0]*R.A[i3]);
    }
    return P;
}
template<class T> HDV kCArray3D<T>* dot(const kCArray4D<T>& L, const kCArray1D<T>& R)
{   kCArray3D<T>* P = new kCArray3D<T>(L.S0,L.E0,L.S1,L.E1,L.S2,L.E2);
    for(int i0=0;i0<L.N0;i0++)
    for(int i1=0;i1<L.N1;i1++)
    for(int i2=0;i2<L.N2;i2++)
    {   P->A[(i0*L.N1+i1)*L.N2+i2]=(T) 0.; // assumes R.N0=L.N3
        for(int i3=0;i3<L.N3;i3++) P->A[(i0*L.N0+i1)*L.N1+i2] += (L.A[((i0*L.N1+i1)*L.N2+i2)*L.N3+i3]*R.A[i3]);
    }
    return P;
}

template<class T> HDV kFArray3D<T>* dot(const kFArray5D<T>& L, const kFArray2D<T>& R)
{   kFArray3D<T>* P = new kFArray3D<T>(L.S0,L.E0,L.S1,L.E1,L.S2,L.E2);
    for(int i0=0;i0<L.N0;i0++)
    for(int i1=0;i1<L.N1;i1++)
    for(int i2=0;i2<L.N2;i2++)
    {   P->A[(i2*L.N1+i1)*L.N0+i0]=(T) 0.; // assumes R.N0=L.N3, R.N1=L.N4
        for(int i4=0;i4<L.N4;i4++)
        for(int i3=0;i3<L.N3;i3++) P->A[(i2*L.N1+i1)*L.N0+i0] += (L.A[(((i4*L.N3+i3)*L.N2+i2)*L.N1+i1)*L.N0+i0]*R.A[i4*L.N3+i3]);
    }
    return P;
}
template<class T> HDV kCArray3D<T>* dot(const kCArray5D<T>& L, const kCArray2D<T>& R)
{   kCArray3D<T>* P = new kCArray3D<T>(L.S0,L.E0,L.S1,L.E1,L.S2,L.E2);
    for(int i0=0;i0<L.N0;i0++)
    for(int i1=0;i1<L.N1;i1++)
    for(int i2=0;i2<L.N2;i2++)
    {   P->A[(i0*L.N1+i1)*L.N2+i2]=(T) 0.; // assumes R.N0=L.N3, R.N1=L.N4
        for(int i3=0;i3<L.N3;i3++)
        for(int i4=0;i4<L.N4;i4++) P->A[(i0*L.N1+i1)*L.N2+i2] += (L.A[(((i0*L.N1+i1)*L.N2+i2)*L.N3+i3)*L.N4+i4]*R.A[i3*L.N4+i4]);
    }
    return P;
}

template<class T> HDV kFArray3D<T>* dot(const kFArray6D<T>& L, const kFArray3D<T>& R)
{   kFArray3D<T>* P = new kFArray3D<T>(L.S0,L.E0,L.S1,L.E1,L.S2,L.E2);
    for(int i0=0;i0<L.N0;i0++)
    for(int i1=0;i1<L.N1;i1++)
    for(int i2=0;i2<L.N2;i2++)
    {   P->A[(i2*L.N1+i1)*L.N0+i0]=(T) 0.; // assumes R.N0=L.N3, R.N1=L.N4, R.N2=L.N5
        for(int i3=0;i3<L.N3;i3++)
        for(int i4=0;i4<L.N4;i4++)
        for(int i5=0;i5<L.N5;i5++) P->A[(i2*L.N2+i1)*L.N1+i0] += (L.A[((((i5*L.N4+i4)*L.N3+i3)*L.N2+i2)*L.N1+i1)*L.N0+i0]*R.A[(i5*L.N4+i4)*L.N3+i3]);
    }
    return P;
}
template<class T> HDV kCArray3D<T>* dot(const kCArray6D<T>& L, const kCArray3D<T>& R)
{   kCArray3D<T>* P = new kCArray3D<T>(L.S0,L.E0,L.S1,L.E1,L.S2,L.E2);
    for(int i0=0;i0<L.N0;i0++)
    for(int i1=0;i1<L.N1;i1++)
    for(int i2=0;i2<L.N2;i2++)
    {   P->A[(i0*L.N1+i1)*L.N2+i2]=(T) 0.; // assumes R.N0=L.N3, R.N1=L.N4, R.N2=L.N5
        for(int i3=0;i3<L.N3;i3++)
        for(int i4=0;i4<L.N4;i4++)
        for(int i5=0;i5<L.N5;i5++) P->A[(i0*L.N1+i1)*L.N2+i2] += (L.A[((((i0*L.N1+i1)*L.N2+i2)*L.N3+i3)*L.N4+i4)*L.N5+i5]*R.A[(i3*L.N4+i4)*L.N5+i5]);
    }
    return P;
}

template<class T> HDV kFArray3D<T>* dot(const kFArray7D<T>& L, const kFArray4D<T>& R)
{   kFArray3D<T>* P = new kFArray3D<T>(L.S0,L.E0,L.S1,L.E1,L.S2,L.E2);
    for(int i0=0;i0<L.N0;i0++)
    for(int i1=0;i1<L.N1;i1++)
    for(int i2=0;i2<L.N2;i2++)
    {   P->A[(i2*L.N1+i1)*L.N0+i0]=(T) 0.; // assumes R.N0=L.N3, R.N1=L.N4, R.N2=L.N5, R.N3=L.N6
        for(int i3=0;i3<L.N3;i3++)
        for(int i4=0;i4<L.N4;i4++)
        for(int i5=0;i5<L.N5;i5++)
        for(int i6=0;i6<L.N6;i6++) P->A[(i2*L.N1+i1)*L.N0+i0] += (L.A[(((((i6*L.N5+i5)*L.N4+i4)*L.N3+i3)*L.N2+i2)*L.N1+i1)*L.N0+i0]*R.A[((i6*L.N5+i5)*L.N4+i4)*L.N3+i3]);
    }
    return P;
}
template<class T> HDV kCArray3D<T>* dot(const kCArray7D<T>& L, const kCArray4D<T>& R)
{   kCArray3D<T>* P = new kCArray3D<T>(L.S0,L.E0,L.S1,L.E1,L.S2,L.E2);
    for(int i0=0;i0<L.N0;i0++)
    for(int i1=0;i1<L.N1;i1++)
    for(int i2=0;i2<L.N2;i2++)
    {   P->A[(i0*L.N1+i1)*L.N2+i2]=(T) 0.; // assumes R.N0=L.N3, R.N1=L.N4, R.N2=L.N5, R.N3=L.N6
        for(int i3=0;i3<L.N3;i3++)
        for(int i4=0;i4<L.N4;i4++)
        for(int i5=0;i5<L.N5;i5++)
        for(int i6=0;i6<L.N6;i6++) P->A[(i0*L.N1+i1)*L.N2+i2] += (L.A[(((((i0*L.N1+i1)*L.N2+i2)*L.N3+i3)*L.N4+i4)*L.N5+i5)*L.N6+i6]*R.A[((i3*L.N4+i4)*L.N5+i5)*L.N6+i6]);
    }
    return P;
}

// dot products returning 4D matrix
template<class T> HDV kFArray4D<T>* dot(const kFArray5D<T>& L, const kFArray1D<T>& R)
{   kFArray4D<T>* P = new kFArray4D<T>(L.S0,L.E0,L.S1,L.E1,L.S2,L.E2,L.S3,L.E3);
    for(int i0=0;i0<L.N0;i0++)
    for(int i1=0;i1<L.N1;i1++)
    for(int i2=0;i2<L.N2;i2++)
    for(int i3=0;i3<L.N3;i3++)
    {   P->A[((i3*L.N2+i2)*L.N1+i1)*L.N0+i0]=(T) 0.; // assumes R.N0=L.N4
        for(int i4=0;i4<L.N4;i4++) P->A[((i3*L.N2+i2)*L.N1+i1)*L.N0+i0] += (L.A[(((i4*L.N3+i3)*L.N2+i2)*L.N1+i1)*L.N0+i0]*R.A[i4]);
    }
    return P;
}
template<class T> HDV kCArray4D<T>* dot(const kCArray5D<T>& L, const kCArray1D<T>& R)
{   kCArray4D<T>* P = new kCArray4D<T>(L.S0,L.E0,L.S1,L.E1,L.S2,L.E2,L.S3,L.E3);
    for(int i0=0;i0<L.N0;i0++)
    for(int i1=0;i1<L.N1;i1++)
    for(int i2=0;i2<L.N2;i2++)
    for(int i3=0;i3<L.N3;i3++)
    {   P->A[((i0*L.N1+i1)*L.N2+i2)*L.N3+i3]=(T) 0.; // assumes R.N0=L.N4
        for(int i4=0;i4<L.N4;i4++) P->A[((i0*L.N1+i1)*L.N2+i2)*L.N3+i3] += (L.A[(((i0*L.N1+i1)*L.N2+i2)*L.N3+i3)*L.N3+i4]*R.A[i4]);
    }
    return P;
}

template<class T> HDV kFArray4D<T>* dot(const kFArray6D<T>& L, const kFArray2D<T>& R)
{   kFArray4D<T>* P = new kFArray4D<T>(L.S0,L.E0,L.S1,L.E1,L.S2,L.E2,L.S3,L.E3);
    for(int i0=0;i0<L.N0;i0++)
    for(int i1=0;i1<L.N1;i1++)
    for(int i2=0;i2<L.N2;i2++)
    for(int i3=0;i3<L.N3;i3++)
    {   P->A[((i3*L.N2+i2)*L.N1+i1)*L.N0+i0]=(T) 0.; // assumes R.N0=L.N4, R.N1=L.N5
        for(int i4=0;i4<L.N4;i4++)
        for(int i5=0;i5<L.N5;i5++) P->A[((i3*L.N2+i2)*L.N1+i1)*L.N0+i0] += (L.A[((((i5*L.N4+i4)*L.N3+i3)*L.N2+i2)*L.N1+i1)*L.N0+i0]*R.A[i5*L.N4+i4]);
    }
    return P;
}
template<class T> HDV kCArray4D<T>* dot(const kCArray6D<T>& L, const kCArray2D<T>& R)
{   kCArray4D<T>* P = new kCArray4D<T>(L.S0,L.E0,L.S1,L.E1,L.S2,L.E2,L.S3,L.E3);
    for(int i0=0;i0<L.N0;i0++)
    for(int i1=0;i1<L.N1;i1++)
    for(int i2=0;i2<L.N2;i2++)
    for(int i3=0;i3<L.N3;i3++)
    {   P->A[((i0*L.N1+i1)*L.N2+i2)*L.N3+i3]=(T) 0.; // assumes R.N0=L.N4, R.N1=L.N5
        for(int i4=0;i4<L.N4;i4++)
        for(int i5=0;i5<L.N5;i5++) P->A[((i0*L.N1+i1)*L.N2+i2)*L.N3+i3] += (L.A[((((i0*L.N1+i1)*L.N2+i2)*L.N3+i3)*L.N3+i4)*L.N5+i5]*R.A[i4*L.N5+i5]);
    }
    return P;
}

template<class T> HDV kFArray4D<T>* dot(const kFArray7D<T>& L, const kFArray3D<T>& R)
{   kFArray4D<T>* P = new kFArray4D<T>(L.S0,L.E0,L.S1,L.E1,L.S2,L.E2,L.S3,L.E3);
    for(int i0=0;i0<L.N0;i0++)
    for(int i1=0;i1<L.N1;i1++)
    for(int i2=0;i2<L.N2;i2++)
    for(int i3=0;i3<L.N3;i3++)
    {   P->A[((i3*L.N2+i2)*L.N1+i1)*L.N0+i0]=(T) 0.; // assumes R.N0=L.N4, R.N1=L.N5, R.N2=L.N6
        for(int i4=0;i4<L.N4;i4++)
        for(int i5=0;i5<L.N5;i5++)
        for(int i6=0;i6<L.N6;i6++) P->A[((i3*L.N2+i2)*L.N1+i1)*L.N0+i0] += (L.A[(((((i6*L.N5+i5)*L.N4+i4)*L.N3+i3)*L.N2+i2)*L.N1+i1)*L.N0+i0]*R.A[(i6*L.N5+i5)*L.N4+i4]);
    }
    return P;
}
template<class T> HDV kCArray4D<T>* dot(const kCArray7D<T>& L, const kCArray3D<T>& R)
{   kFArray4D<T>* P = new kFArray4D<T>(L.S0,L.E0,L.S1,L.E1,L.S2,L.E2,L.S3,L.E3);
    for(int i0=0;i0<L.N0;i0++)
    for(int i1=0;i1<L.N1;i1++)
    for(int i2=0;i2<L.N2;i2++)
    for(int i3=0;i3<L.N3;i3++)
    {   P->A[((i0*L.N1+i1)*L.N2+i2)*L.N3+i3]=(T) 0.; // assumes R.N0=L.N4, R.N1=L.N5
        for(int i4=0;i4<L.N4;i4++)
        for(int i5=0;i5<L.N5;i5++)
        for(int i6=0;i6<L.N6;i6++) P->A[((i0*L.N1+i1)*L.N2+i2)*L.N3+i3] += (L.A[(((((i0*L.N1+i1)*L.N2+i2)*L.N3+i3)*L.N3+i4)*L.N5+i5)*L.N6+i6]*R.A[(i4*L.N5+i5)*L.N6+i6]);
    }
    return P;
}

// dot products returning 5D matrix
template<class T> HDV kFArray5D<T>* dot(const kFArray6D<T>& L, const kFArray1D<T>& R)
{   kFArray5D<T>* P = new kFArray5D<T>(L.S0,L.E0,L.S1,L.E1,L.S2,L.E2,L.S3,L.E3,L.S4,L.E4);
    for(int i0=0;i0<L.N0;i0++)
    for(int i1=0;i1<L.N1;i1++)
    for(int i2=0;i2<L.N2;i2++)
    for(int i3=0;i3<L.N3;i3++)
    for(int i4=0;i4<L.N4;i4++)
    {   P->A[(((i4*L.N3+i3)*L.N2+i2)*L.N1+i1)*L.N0+i0]=(T) 0.; // assumes R.N0=L.N5
        for(int i5=0;i5<L.N5;i5++) P->A[(((i4*L.N3+i3)*L.N2+i2)*L.N1+i1)*L.N0+i0] += (L.A[((((i5*L.N4+i4)*L.N3+i3)*L.N2+i2)*L.N1+i1)*L.N0+i0]*R.A[i5]);
    }
    return P;
}
template<class T> HDV kCArray5D<T>* dot(const kCArray6D<T>& L, const kCArray1D<T>& R)
{   kCArray5D<T>* P = new kCArray5D<T>(L.S0,L.E0,L.S1,L.E1,L.S2,L.E2,L.S3,L.E3,L.S4,L.E4);
    for(int i0=0;i0<L.N0;i0++)
    for(int i1=0;i1<L.N1;i1++)
    for(int i2=0;i2<L.N2;i2++)
    for(int i3=0;i3<L.N3;i3++)
    for(int i4=0;i4<L.N4;i4++)
    {   P->A[(((i0*L.N1+i1)*L.N2+i2)*L.N3+i3)*L.N4+i4]=(T) 0.; // assumes R.N0=L.N5
        for(int i5=0;i5<L.N5;i5++) P->A[(((i0*L.N1+i1)*L.N2+i2)*L.N3+i3)*L.N4+i4] += (L.A[((((i0*L.N1+i1)*L.N2+i2)*L.N3+i3)*L.N3+i4)*L.N5+i5]*R.A[i5]);
    }
    return P;
}

template<class T> HDV kFArray5D<T>* dot(const kFArray7D<T>& L, const kFArray2D<T>& R)
{   kFArray5D<T>* P = new kFArray5D<T>(L.S0,L.E0,L.S1,L.E1,L.S2,L.E2,L.S3,L.E3,L.S4,L.E4);
    for(int i0=0;i0<L.N0;i0++)
    for(int i1=0;i1<L.N1;i1++)
    for(int i2=0;i2<L.N2;i2++)
    for(int i3=0;i3<L.N3;i3++)
    for(int i4=0;i4<L.N4;i4++)
    {   P->A[(((i4*L.N3+i3)*L.N2+i2)*L.N1+i1)*L.N0+i0]=(T) 0.; // assumes R.N0=L.N5, R.N1=L.N6
        for(int i5=0;i5<L.N5;i5++)
        for(int i6=0;i6<L.N6;i6++) P->A[(((i4*L.N3+i3)*L.N2+i2)*L.N1+i1)*L.N0+i0] += (L.A[(((((i6*L.N5+i5)*L.N4+i4)*L.N3+i3)*L.N2+i2)*L.N1+i1)*L.N0+i0]*R.A[i6*L.N5+i5]);
    }
    return P;
}
template<class T> HDV kCArray5D<T>* dot(const kCArray7D<T>& L, const kCArray2D<T>& R)
{   kFArray5D<T>* P = new kFArray5D<T>(L.S0,L.E0,L.S1,L.E1,L.S2,L.E2,L.S3,L.E3,L.S4,L.E4);
    for(int i0=0;i0<L.N0;i0++)
    for(int i1=0;i1<L.N1;i1++)
    for(int i2=0;i2<L.N2;i2++)
    for(int i3=0;i3<L.N3;i3++)
    for(int i4=0;i4<L.N4;i4++)
    {   P->A[(((i0*L.N1+i1)*L.N2+i2)*L.N3+i3)*L.N4+i4]=(T) 0.; // assumes R.N0=L.N5, R.N1=L.N6
        for(int i5=0;i5<L.N5;i5++)
        for(int i6=0;i6<L.N6;i6++) P->A[(((i0*L.N1+i1)*L.N2+i2)*L.N3+i3)*L.N4+i4] += (L.A[(((((i0*L.N1+i1)*L.N2+i2)*L.N3+i3)*L.N3+i4)*L.N5+i5)*L.N6+i6]*R.A[i5*L.N6+i6]);
    }
    return P;
}

// dot products returning 6D matrix
template<class T> HDV kFArray6D<T>& dot(const kFArray7D<T>& L, const kFArray1D<T>& R)
{   kFArray6D<T>* P = new kFArray6D<T>(L.S0,L.E0,L.S1,L.E1,L.S2,L.E2,L.S3,L.E3,L.S4,L.E4,L.S5,L.E5);
    for(int i0=0;i0<L.N0;i0++)
    for(int i1=0;i1<L.N1;i1++)
    for(int i2=0;i2<L.N2;i2++)
    for(int i3=0;i3<L.N3;i3++)
    for(int i4=0;i4<L.N4;i4++)
    for(int i5=0;i5<L.N5;i5++)
    {   P->A[((((i5*L.N4+i4)*L.N3+i3)*L.N2+i2)*L.N1+i1)*L.N0+i0]=(T) 0.; // assumes R.N0=L.N6
        for(int i6=0;i6<L.N6;i6++) P->A[((((i5*L.N4+i4)*L.N3+i3)*L.N2+i2)*L.N1+i1)*L.N0+i0] += (L.A[(((((i6*L.N5+i5)*L.N4+i4)*L.N3+i3)*L.N2+i2)*L.N1+i1)*L.N0+i0]*R.A[i6]);
    }
    return P;
}
template<class T> HDV kCArray6D<T>* dot(const kCArray7D<T>& L, const kCArray1D<T>& R)
{   kCArray6D<T>* P = new kCArray6D<T>(L.S0,L.E0,L.S1,L.E1,L.S2,L.E2,L.S3,L.E3,L.S4,L.E4,L.S5,L.E5);
    for(int i0=0;i0<L.N0;i0++)
    for(int i1=0;i1<L.N1;i1++)
    for(int i2=0;i2<L.N2;i2++)
    for(int i3=0;i3<L.N3;i3++)
    for(int i4=0;i4<L.N4;i4++)
    for(int i5=0;i5<L.N5;i5++)
    {   P->A[((((i0*L.N1+i1)*L.N2+i2)*L.N3+i3)*L.N4+i4)*L.N5+i5]=(T) 0.; // assumes R.N0=L.N6
        for(int i6=0;i6<L.N6;i6++) P->A[((((i0*L.N1+i1)*L.N2+i2)*L.N3+i3)*L.N4+i4)*L.N5+i5] += (L.A[(((((i0*L.N1+i1)*L.N2+i2)*L.N3+i3)*L.N3+i4)*L.N5+i5)*L.N6+i6]*R.A[i6]);
    }
    return P;
}


// Linear Algebra

template <class T> HDV void implicitSolve(kFArray1D<T>* pS, kFArray2D<T>* pM, kFArray1D<T>* pR)
{
    //  IMPORTANT NOTE: A[][] and rhs[] are MODIFIED !!
    int np,ns,ms;
    int N = pR->N0;
    if((N != pM->N1) || (N != pM->N0) || (N != pS->N0))
    {   printf("Kernel::implicitSolve::error: dimensions do not match\n");
        return;
    }
    int Nmax = N-1;

    T errP;
    T* jac = pM->A;
    T* sol = pS->A;
    T* rhs = pR->A;

    for(ns=0;ns<=Nmax;ns++) sol[ns]=rhs[ns];
    for(np=0;np< Nmax;np++)
    {
        for(ms=np+1;ms<=Nmax;ms++)
        {   errP = jac[np*N+ms]/jac[np*N+np];
            for(ns=np+1;ns<=Nmax;ns++) jac[ns*N+ms] -= (errP*jac[ns*N+np]);
            jac[np*N+ms] = 0.;
            sol[ms] -= errP*sol[np];
        }
    }
    sol[Nmax] = sol[Nmax]/jac[Nmax*N+Nmax];
    for(np=Nmax-1;np>=0;np--)
    {
        errP=0.;
        for(ns=np+1;ns<=Nmax;ns++) errP += (jac[ns*N+np]*sol[ns]);
        sol[np] = (sol[np]-errP)/jac[np*N+np] ;
    }
} // end Point-Implicit Solver


template <class T> HDV void implicitSolve(kCArray1D<T>* pS, kCArray2D<T>* pM, kCArray1D<T>* pR)
{
    //  IMPORTANT NOTE: A[][] and rhs[] are MODIFIED !!
    int np,ns,ms;
    int N = pR->N0;
    if((N != pM->N1) || (N != pM->N0) || (N != pS->N0))
    {   printf("implicitSolve::error: dimensions do not match\n");
        return;
    }
    int Nmax = N-1;
    T errP;
    T* jac = pM->A;
    T* sol = pS->A;
    T* rhs = pR->A;

    for(ns=0;ns<=Nmax;ns++) sol[ns]=rhs[ns];
    for(np=0;np< Nmax;np++)
    {
        for(ms=np+1;ms<=Nmax;ms++)
        {   errP = jac[ms*N+np]/jac[np*N+np];
            for(ns=np+1;ns<=Nmax;ns++) jac[ms*N+ns] -= (errP*jac[np*N+ns]);
            jac[ms*N+np] = 0.;
            sol[ms] -= errP*sol[np];
        }
    }
    sol[Nmax] = sol[Nmax]/jac[Nmax*N+Nmax];
    for(np=Nmax-1;np>=0;np--)
    {
        errP=0.;
        for(ns=np+1;ns<=Nmax;ns++) errP += (jac[np*N+ns]*sol[ns]);
        sol[np] = (sol[np]-errP)/jac[np*N+np] ;
    }
} // end Point-Implicit Solver




#endif	/* KMATRIX_H */

