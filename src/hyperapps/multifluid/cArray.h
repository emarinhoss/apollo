/*
 * cArray.h
 *
 *  Created on: Oct 9, 2016
 *      Author: hle
 */

#ifndef CARRAY_H_
#define CARRAY_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cmath>
#include <time.h>
#include <vector>
#include "constant.h"

// C-array utilities
using namespace std;
// 1D
template<class T>  T* cArray(int N0)
{
    T* A = new T[N0];
    for(int i=0;i<N0;i++) A[i] = (T)0;
    return A;
}
template<class T> void del_cArray(int N0, T* A)
{
    if(A)
    {
        delete[] A;
    }
}
template<class T> void set_cArray(int N0, T* A, T X)
{
	for(int i=0;i<N0;i++) A[i] = (T)X;
}
template<class T> void cpy_cArray(int N0, T* A, T* B)
{
	for(int i=0;i<N0;i++) A[i] = B[i];
}

// 2D
template<class T>  T** cArray(int N0, int N1)
{
    T** A = new T*[N0];
    for(int i=0;i<N0;i++)
    {
        A[i] = new T[N1];
        for(int j=0;j<N1;j++) A[i][j] = (T)0;
    }
    return A;
}
template<class T> void del_cArray(int N0, int N1, T** A)
{
    if(A)
    {
        for(int i=0;i<N0;i++) if(A[i]){delete[] A[i];A[i]=NULL;}
        delete[] A;A=NULL;
    }
}
template<class T> void set_cArray(int N0, int N1, T** A, T X)
{
    for(int i=0;i<N0;i++)
    for(int j=0;j<N1;j++)
    	A[i][j] = (T)X;
}
template<class T> void cpy_cArray(int N0, int N1, T** A, T** B)
{
	for(int i=0;i<N0;i++)
	for(int j=0;j<N1;j++)
		A[i][j] = B[i][j];

}


// 3D
template<class T>  T*** cArray(int N0, int N1, int N2)
{
    T*** A = new T**[N0];
    for(int i=0;i<N0;i++)
    {
        A[i] = new T*[N1];
        for(int j=0;j<N1;j++)
        {
            A[i][j] = new T[N2];
            for(int k=0;k<N2;k++) A[i][j][k] = (T)0;
        }
    }
    return A;
}
template<class T> void del_cArray(int N0, int N1, int N2, T*** A)
{
    if(A)
    {
        for(int i=0;i<N0;i++)
        {
            for(int j=0;j<N1;j++) delete[] A[i][j];
            delete[] A[i];
        }
        delete[] A;
    }
}
template<class T> void set_cArray(int N0, int N1, int N2, T*** A, T X)
{
    for(int i=0;i<N0;i++)
    for(int j=0;j<N1;j++)
	for(int k=0;k<N2;k++)
    	A[i][j][k] = (T)X;
}
template<class T> void cpy_cArray(int N0, int N1, int N2, T*** A, T*** B)
{
	for(int i=0;i<N0;i++)
	for(int j=0;j<N1;j++)
	for(int k=0;k<N2;k++)
		A[i][j][k] = B[i][j][k];

}


// 4D
template<class T>  T**** cArray(int N0, int N1, int N2, int N3)
{
    T**** A = new T***[N0];
    for(int i=0;i<N0;i++)
    {
        A[i] = new T**[N1];
        for(int j=0;j<N1;j++)
        {
            A[i][j] = new T*[N2];
            for(int k=0;k<N2;k++)
            {
                A[i][j][k] = new T[N3];
                for(int m=0;m<N3;m++) A[i][j][k][m] = (T)0;
            }
        }
    }
    return A;
}
template<class T> void del_cArray(int N0, int N1, int N2, int N3, T**** A)
{
    if(A)
    {
        for(int i=0;i<N0;i++)
        {
            for(int j=0;j<N1;j++)
            {
                for(int k=0;k<N2;k++) delete[] A[i][j][k];
                delete[] A[i][j];
            }
            delete[] A[i];
        }
        delete[] A;
    }
}
template<class T> void set_cArray(int N0, int N1, int N2, int N3, T**** A, T X)
{
    for(int i=0;i<N0;i++)
    for(int j=0;j<N1;j++)
	for(int k=0;k<N2;k++)
	for(int l=0;l<N3;l++)
    	A[i][j][k][l] = (T)X;
}
template<class T> void cpy_cArray(int N0, int N1, int N2, int N3, T**** A, T**** B)
{
	for(int i=0;i<N0;i++)
	for(int j=0;j<N1;j++)
	for(int k=0;k<N2;k++)
	for(int l=0;l<N3;l++)
		A[i][j][k][l] = B[i][j][k][l];

}


// 5D
template<class T>  T***** cArray(int N0, int N1, int N2, int N3, int N4)
{
    T***** A = new T****[N0];
    for(int i=0;i<N0;i++)
    {
        A[i] = new T***[N1];
        for(int j=0;j<N1;j++)
        {
            A[i][j] = new T**[N2];
            for(int k=0;k<N2;k++)
            {
                A[i][j][k] = new T*[N3];
                for(int m=0;m<N3;m++)
                {
                    A[i][j][k][m] = new T[N4];
                    for(int n=0;n<N4;n++) A[i][j][k][m][n] = (T)0;
                }
            }
        }
    }
    return A;
}
template<class T> void del_cArray(int N0, int N1, int N2, int N3, int N4, T***** A)
{
    if(A)
    {
        for(int i=0;i<N0;i++)
        {
            for(int j=0;j<N1;j++)
            {
                for(int k=0;k<N2;k++)
                {
                    for(int m=0;m<N3;m++) delete[] A[i][j][k][m];
                    delete[] A[i][j][k];
                }
                delete[] A[i][j];
            }
            delete[] A[i];
        }
        delete[] A;
    }
}
template<class T> void set_cArray(int N0, int N1, int N2, int N3, int N4, T***** A, T X)
{
    for(int i=0;i<N0;i++)
    for(int j=0;j<N1;j++)
	for(int k=0;k<N2;k++)
	for(int l=0;l<N3;l++)
	for(int n=0;n<N4;n++)
    	A[i][j][k][l][n] = (T)X;
}
template<class T> void cpy_cArray(int N0, int N1, int N2, int N3, int N4, T***** A, T***** B)
{
	for(int i=0;i<N0;i++)
	for(int j=0;j<N1;j++)
	for(int k=0;k<N2;k++)
	for(int l=0;l<N3;l++)
	for(int n=0;n<N4;n++)
		A[i][j][k][l][n] = B[i][j][k][l][n];

}

// 6D
template<class T>  T****** cArray(int N0, int N1, int N2, int N3, int N4, int N5)
{
    T****** A = new T*****[N0];
    for(int i=0;i<N0;i++)
    {
        A[i] = new T****[N1];
        for(int j=0;j<N1;j++)
        {
            A[i][j] = new T***[N2];
            for(int k=0;k<N2;k++)
            {
                A[i][j][k] = new T**[N3];
                for(int m=0;m<N3;m++)
                {
                    A[i][j][k][m] = new T*[N4];
                    for(int n=0;n<N4;n++)
                    {
                        A[i][j][k][m][n] = new T[N5];
                        for(int p=0;p<N5;p++) A[i][j][k][m][n][p] = (T)0;
                    }
                }
            }
        }
    }
    return A;
}
template<class T> void del_cArray(int N0, int N1, int N2, int N3, int N4, int N5, T****** A)
{
    if(A)
    {
        for(int i=0;i<N0;i++)
        {
            for(int j=0;j<N1;j++)
            {
                for(int k=0;k<N2;k++)
                {
                    for(int m=0;m<N3;m++)
                    {
                        for(int n=0;n<N4;n++) delete[] A[i][j][k][m][n];
                        delete[] A[i][j][k][m];
                    }
                    delete[] A[i][j][k];
                }
                delete[] A[i][j];
            }
            delete[] A[i];
        }
        delete[] A;
    }
}
template<class T> void set_cArray(int N0, int N1, int N2, int N3, int N4, int N5, T****** A, T X)
{
    for(int i=0;i<N0;i++)
    for(int j=0;j<N1;j++)
	for(int k=0;k<N2;k++)
	for(int l=0;l<N3;l++)
	for(int n=0;n<N4;n++)
	for(int o=0;o<N5;o++)
    	A[i][j][k][l][n][o] = (T)X;
}
template<class T> void cpy_cArray(int N0, int N1, int N2, int N3, int N4, int N5, T****** A, T****** B)
{
	for(int i=0;i<N0;i++)
	for(int j=0;j<N1;j++)
	for(int k=0;k<N2;k++)
	for(int l=0;l<N3;l++)
	for(int n=0;n<N4;n++)
	for(int o=0;o<N5;o++)
		A[i][j][k][l][n][o] = B[i][j][k][l][n][o];

}

// 7D
template<class T>  T******* cArray(int N0, int N1, int N2, int N3, int N4, int N5, int N6)
{
    T******* A = new T******[N0];
    for(int i=0;i<N0;i++)
    {
        A[i] = new T*****[N1];
        for(int j=0;j<N1;j++)
        {
            A[i][j] = new T****[N2];
            for(int k=0;k<N2;k++)
            {
                A[i][j][k] = new T***[N3];
                for(int m=0;m<N3;m++)
                {
                    A[i][j][k][m] = new T**[N4];
                    for(int n=0;n<N4;n++)
                    {
                        A[i][j][k][m][n] = new T*[N5];
                        for(int p=0;p<N5;p++)
                        {
                            A[i][j][k][m][n][p] = new T[N6];
                            for(int q=0;q<N6;q++) A[i][j][k][m][n][p][q] = (T)0;
                        }
                    }
                }
            }
        }
    }
    return A;
}
template<class T> void del_cArray(int N0, int N1, int N2, int N3, int N4, int N5, int N6, T******* A)
{
    if(A)
    {
        for(int i=0;i<N0;i++)
        {
            for(int j=0;j<N1;j++)
            {
                for(int k=0;k<N2;k++)
                {
                    for(int m=0;m<N3;m++)
                    {
                        for(int n=0;n<N4;n++)
                        {
                            for(int p=0;p<N5;p++) delete[] A[i][j][k][m][n][p];
                            delete[] A[i][j][k][m][n];
                        }
                        delete[] A[i][j][k][m];
                    }
                    delete[] A[i][j][k];
                }
                delete[] A[i][j];
            }
            delete[] A[i];
        }
        delete[] A;
    }
}
template<class T> void set_cArray(int N0, int N1, int N2, int N3, int N4, int N5, int N6, T******* A, T X)
{
    for(int i=0;i<N0;i++)
    for(int j=0;j<N1;j++)
	for(int k=0;k<N2;k++)
	for(int l=0;l<N3;l++)
	for(int n=0;n<N4;n++)
	for(int o=0;o<N5;o++)
	for(int p=0;p<N6;p++)
    	A[i][j][k][l][n][o][p] = (T)X;
}
template<class T> void cpy_cArray(int N0, int N1, int N2, int N3, int N4, int N5, int N6, T******* A, T******* B)
{
	for(int i=0;i<N0;i++)
	for(int j=0;j<N1;j++)
	for(int k=0;k<N2;k++)
	for(int l=0;l<N3;l++)
	for(int n=0;n<N4;n++)
	for(int o=0;o<N5;o++)
	for(int p=0;p<N6;p++)
		A[i][j][k][l][n][o][p] = B[i][j][k][l][n][o][p];

}

template<class T> T* linspace(double _low, double _upp, int N)
{
	T _d = (_upp-_low)/(T)(N-1);
	T* p = new T[N];
	for(int i=0;i<N;i++) p[i] = _low+(double)i*_d;
	return p;
}

template<class T> T* logspace(double _low, double _upp, int N)
{
	T _d = (log(_upp)-log(_low))/(T)(N-1);
	T* p = new T[N];
	for(int i=0;i<N;i++) p[i] = exp(log(_low)+(double)i*_d);
	return p;
}





#endif /* CARRAY_H_ */
