#ifndef WXPNODALDGFUNCTIONS_H
#define WXPNODALDGFUNCTIONS_H

#include <wxlogger.h>
#include <wxlogstream.h>
#include "wxmath.h"

template<typename REAL>
void
dataN01(int N, REAL *pr, REAL *ps, REAL *pDr, REAL *pDs, REAL *pDrw, REAL *pDsw, REAL *pVand, REAL *pIVand, REAL *pLIFT, int *pFmask)
{
    #include "dataN01.h"
    int NpE = (N+1)*(N+2)/2;
    int NpF = (N+1);
    int NfE = 3;

    /* load r, s */
    for(int n=0;n<NpE;++n){
      pr[n] = p_r[n];
      ps[n] = p_s[n];
    }

    int sk = 0;
    /* load Dr, Ds, Drw, Drs, Vandermonde and Inverse Vandermonde */
    for(int n=0;n<NpE;++n){
      for(int m=0;m<NpE;++m){
        pVand[sk] = p_Vand[n][m];
        pIVand[sk] = p_IVand[n][m];
        pDr[sk] = p_Dr[n][m];
        pDs[sk] = p_Ds[n][m];
        pDrw[sk] = p_Drw[n][m];
        pDsw[sk++] = p_Dsw[n][m];
      }
    }

    sk = 0;
    /* load LIFT */
    for(int n=0;n<NpE;++n){
      for(int m=0;m<NpF*NfE;++m){
        pLIFT[sk++] = p_LIFT[n][m];
      }
    }

    sk = 0;
    for(int n=0;n<NfE;++n){
      for(int m=0;m<NpF;++m){
        pFmask[sk++] = p_Fmask[n][m];
      }
    }
}



template<typename REAL>
void
nodalDGfunctions(unsigned N, REAL *p_r, REAL *p_s, REAL *p_Dr, REAL *p_Ds, REAL *p_Drw, REAL *p_Dsw, REAL *p_Vand, REAL *p_IVand, REAL *p_LIFT, int *p_Fmask)
{
    switch (N) {
    case 1:
        dataN01(N, p_r, p_s, p_Dr, p_Ds , p_Drw, p_Dsw, p_Vand, p_IVand, p_LIFT, p_Fmask);
        break;
//    case 2:
//        dataN02(p_r, p_s, p_Dr, p_Ds, p_LIFT, p_Fmask);
//        break;
//    case 3:
//        dataN03(p_r, p_s, p_Dr, p_Ds, p_LIFT, p_Fmask);
//        break;
//    case 4:
//        dataN04(p_r, p_s, p_Dr, p_Ds, p_LIFT, p_Fmask);
//        break;
//    case 5:
//        dataN05(p_r, p_s, p_Dr, p_Ds, p_LIFT, p_Fmask);
//        break;
//    case 6:
//        dataN06(p_r, p_s, p_Dr, p_Ds, p_LIFT, p_Fmask);
//        break;
//    case 7:
//        dataN07(p_r, p_s, p_Dr, p_Ds, p_LIFT, p_Fmask);
//        break;
//    case 8:
//        dataN08(p_r, p_s, p_Dr, p_Ds, p_LIFT, p_Fmask);
//        break;
//    case 9:
//        dataN09(p_r, p_s, p_Dr, p_Ds, p_LIFT, p_Fmask);
//        break;
//    case 10:
//        dataN10(p_r, p_s, p_Dr, p_Ds, p_LIFT, p_Fmask);
//        break;
//    case 11:
//        dataN11(p_r, p_s, p_Dr, p_Ds, p_LIFT, p_Fmask);
//        break;
//    case 12:
//        dataN12(p_r, p_s, p_Dr, p_Ds, p_LIFT, p_Fmask);
//        break;
//    case 13:
//        dataN13(p_r, p_s, p_Dr, p_Ds, p_LIFT, p_Fmask);
//        break;
//    case 14:
//        dataN14(p_r, p_s, p_Dr, p_Ds, p_LIFT, p_Fmask);
//        break;
//    case 15:
//        dataN15(p_r, p_s, p_Dr, p_Ds, p_LIFT, p_Fmask);
//        break;
//    case 16:
//        dataN16(p_r, p_s, p_Dr, p_Ds, p_LIFT, p_Fmask);
//        break;
    default:
        break;
    }
}

template<typename REAL>
void
polyms_Ord_1(int po, REAL *pr, REAL *ps, int *pFmask)
{
    REAL p_r[3] = {-1. , 1 ,-1 };
    REAL p_s[3] = {-1 ,-1 , 1 };
    int p_Fmask[6] = {0,1,1,2,2,0};

    for(unsigned k=0; k<3; k++)
    {
        pr[k] = p_r[k];
        ps[k] = p_s[k];
    }

    for(unsigned k=0; k<6; k++)
        pFmask[k] = p_Fmask[k];

}

template<typename REAL>
void
polyms_Ord_2(int po, REAL *pr, REAL *ps, int *pFmask)
{
    REAL p_r[6] = {                -1 ,                 0 ,                 1 ,                -1 ,                 0 ,                -1 };
    REAL p_s[6] = {               -1 ,                -1 ,                -1 ,                 0 ,                 0 ,                 1 };
    int p_Fmask[9] = {0 , 1 , 2 , 2 , 4 , 5 ,0 , 3 , 5 };

    for(unsigned k=0; k<(po+1)*(po+2)/2; k++)
    {
        pr[k] = p_r[k];
        ps[k] = p_s[k];
    }

    for(unsigned k=0; k<3*(po+1); k++)
        pFmask[k] = p_Fmask[k];

}

template<typename REAL>
void
nodalNaturalCoordinates(unsigned N, REAL *pr, REAL *ps, int *pFmask)
{
    switch (N) {
    case 1: polyms_Ord_1(N,pr,ps,pFmask); break;
    case 2: polyms_Ord_2(N,pr,ps,pFmask); break;
    default:
        WxLogger *log = WxLogger::get("apollo-root.console");
        WxLogStream infStrm = log->getInfoStream();
        infStrm << "** The polynomial order selected " << N << " is not valid. Max N should be 2." << std::endl;
        break;
    }
}

#endif // WXPNODALDGFUNCTIONS_H
