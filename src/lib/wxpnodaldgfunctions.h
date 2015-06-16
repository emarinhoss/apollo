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
    int p_Fmask[9] = {0 , 1 , 2 , 2 , 4 , 5 ,5 , 3 , 0 };

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
polyms_Ord_3(int po, REAL *pr, REAL *ps, int *pFmask)
{
    REAL p_r[10] = {                -1 , -0.447213595499958 , 0.447213595499958 ,                 1 ,                -1 , -0.333333333333333 , 0.447213595499958 ,                -1 , -0.447213595499958 ,                -1 };
    REAL p_s[10] = {               -1 ,                -1 ,                -1 ,                -1 , -0.447213595499958 , -0.333333333333333 , -0.447213595499958 , 0.447213595499958 , 0.447213595499958 ,                 1 };
    int p_Fmask[12] = {0 , 1 , 2 , 3, 3 , 6 , 8 , 9, 9 , 7, 4, 0};

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
polyms_Ord_4(int po, REAL *pr, REAL *ps, int *pFmask)
{
    REAL p_r[15] = {                -1 , -0.654653670707977 ,                 0 , 0.654653670707977 ,                 1 ,                -1 , -0.551583507555306 , 0.103167015110611 , 0.654653670707977 ,                -1 , -0.551583507555306 ,                 0 ,                -1 , -0.654653670707977 ,                -1 };
    REAL p_s[15] = {               -1 ,                -1 ,                -1 ,                -1 ,                -1 , -0.654653670707977 , -0.551583507555305 , -0.551583507555306 , -0.654653670707977 ,                 0 , 0.103167015110611 ,                 0 , 0.654653670707977 , 0.654653670707977 ,                 1 };
    int p_Fmask[15] = {0 , 1 , 2 , 3 , 4, 4 , 8 , 11 , 13 , 14, 14, 12, 9, 5, 0};

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
polyms_Ord_5(int po, REAL *pr, REAL *ps, int *pFmask)
{
    REAL p_r[21] = {                -1 , -0.765055323929465 , -0.285231516480645 , 0.285231516480645 , 0.765055323929465 ,                 1 ,                -1 , -0.684472514501909 , -0.171245477332074 , 0.368945029003817 , 0.765055323929465 ,                -1 , -0.657509045335851 , -0.171245477332074 , 0.285231516480645 ,                -1 , -0.684472514501909 , -0.285231516480645 ,                -1 , -0.765055323929465 ,                -1 };
    REAL p_s[21] = {               -1 ,                -1 ,                -1 ,                -1 ,                -1 ,                -1 , -0.765055323929465 , -0.684472514501909 , -0.657509045335851 , -0.684472514501909 , -0.765055323929465 , -0.285231516480645 , -0.171245477332074 , -0.171245477332074 , -0.285231516480645 , 0.285231516480645 , 0.368945029003817 , 0.285231516480645 , 0.765055323929465 , 0.765055323929465 ,                 1 };
    int p_Fmask[18] = {0 , 1 , 2 , 3 , 4 , 5, 5 , 10 , 14 , 17 , 19 , 20 , 20 , 18 , 15 , 11 , 6 , 0};

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
polyms_Ord_6(int po, REAL *pr, REAL *ps, int *pFmask)
{
    REAL p_r[28] = {                -1 , -0.830223896278567 , -0.468848793470714 ,                 0 , 0.468848793470714 , 0.830223896278567 ,                 1 ,                -1 , -0.769516496193442 , -0.364769256536689 ,  0.10886465257101 , 0.539032992386883 , 0.830223896278567 ,                -1 , -0.744095396034322 , -0.333333333333333 ,  0.10886465257101 , 0.468848793470714 ,                -1 , -0.744095396034321 , -0.364769256536689 ,                 0 ,                -1 , -0.769516496193441 , -0.468848793470714 ,                -1 , -0.830223896278567 ,                -1 };
    REAL p_s[28] = {               -1 ,                -1 ,                -1 ,                -1 ,                -1 ,                -1 ,                -1 , -0.830223896278567 , -0.769516496193442 , -0.744095396034321 , -0.744095396034321 , -0.769516496193441 , -0.830223896278567 , -0.468848793470714 , -0.364769256536689 , -0.333333333333333 , -0.364769256536689 , -0.468848793470714 ,                 0 ,  0.10886465257101 ,  0.10886465257101 ,                 0 , 0.468848793470714 , 0.539032992386883 , 0.468848793470714 , 0.830223896278567 , 0.830223896278567 ,                 1 };
    int p_Fmask[21] = {0 , 1 , 2 , 3 , 4 , 5 , 6 , 6 , 12 , 17 , 21 , 24 , 26 , 27 , 27 , 25 , 22, 18 , 13 , 7 , 0};

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
polyms_Ord_7(int po, REAL *pr, REAL *ps, int *pFmask)
{
    REAL p_r[36] = {                -1 , -0.871740148509607 , -0.591700181433142 , -0.209299217902479 , 0.209299217902479 , 0.591700181433142 , 0.871740148509607 ,                 1 ,                -1 , -0.823996159427635 , -0.503284033260976 , -0.102529737209732 , 0.304713444152606 , 0.647992318855269 , 0.871740148509606 ,                -1 , -0.80142941089163 , -0.465073294429312 , -0.0698534111413768 , 0.304713444152605 , 0.591700181433142 ,                -1 , -0.794940525580537 , -0.465073294429312 , -0.102529737209732 , 0.209299217902479 ,                -1 , -0.80142941089163 , -0.503284033260976 , -0.209299217902479 ,                -1 , -0.823996159427634 , -0.591700181433142 ,                -1 , -0.871740148509606 ,                -1 };
    REAL p_s[36] = {               -1 ,                -1 ,                -1 ,                -1 ,                -1 ,                -1 ,                -1 ,                -1 , -0.871740148509607 , -0.823996159427635 , -0.80142941089163 , -0.794940525580537 , -0.80142941089163 , -0.823996159427635 , -0.871740148509607 , -0.591700181433142 , -0.503284033260976 , -0.465073294429312 , -0.465073294429312 , -0.503284033260976 , -0.591700181433142 , -0.209299217902479 , -0.102529737209732 , -0.0698534111413768 , -0.102529737209732 , -0.209299217902479 , 0.209299217902479 , 0.304713444152606 , 0.304713444152606 , 0.209299217902479 , 0.591700181433142 , 0.647992318855269 , 0.591700181433142 , 0.871740148509607 , 0.871740148509606 ,                 1 };
    int p_Fmask[24] = {0 , 1 , 2 , 3 , 4 , 5 , 6 , 7 , 7 , 14 , 20 , 25 , 29 , 32 , 34 , 35 , 35 , 33 , 30 , 26 , 21 , 15 , 8 , 0};

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
polyms_Ord_8(int po, REAL *pr, REAL *ps, int *pFmask)
{
    REAL p_r[45] = {                -1 , -0.89975799541146 , -0.677186279510738 , -0.363117463826178 ,                 0 , 0.363117463826178 , 0.677186279510738 ,  0.89975799541146 ,                 1 ,                -1 , -0.861819860276539 , -0.602543502716999 , -0.264723799928576 , 0.101891561713898 , 0.446427110295541 , 0.723639720553079 ,  0.89975799541146 ,                -1 , -0.843883607578542 , -0.566537963017939 , -0.221988381746101 , 0.133075926035877 , 0.446427110295541 , 0.677186279510738 ,                -1 , -0.837167761785322 , -0.556023236507798 , -0.221988381746101 , 0.101891561713898 , 0.363117463826178 ,                -1 , -0.837167761785322 , -0.566537963017939 , -0.264723799928576 ,                 0 ,                -1 , -0.843883607578542 , -0.602543502716999 , -0.363117463826178 ,                -1 , -0.861819860276539 , -0.677186279510738 ,                -1 , -0.89975799541146 ,                -1 };
    REAL p_s[45] = {               -1 ,                -1 ,                -1 ,                -1 ,                -1 ,                -1 ,                -1 ,                -1 ,                -1 , -0.89975799541146 , -0.861819860276539 , -0.843883607578542 , -0.837167761785322 , -0.837167761785322 , -0.843883607578542 , -0.861819860276539 , -0.89975799541146 , -0.677186279510738 , -0.602543502716999 , -0.566537963017939 , -0.556023236507798 , -0.566537963017939 , -0.602543502716999 , -0.677186279510738 , -0.363117463826178 , -0.264723799928576 , -0.221988381746101 , -0.221988381746101 , -0.264723799928575 , -0.363117463826178 ,                 0 , 0.101891561713898 , 0.133075926035877 , 0.101891561713898 ,                 0 , 0.363117463826178 , 0.446427110295541 , 0.446427110295541 , 0.363117463826178 , 0.677186279510738 , 0.723639720553079 , 0.677186279510738 ,  0.89975799541146 ,  0.89975799541146 ,                 1 };
    int p_Fmask[27] = {0 , 1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 , 8 , 16 , 23 , 29 , 34 , 38 , 41 , 43 , 44 , 44 , 42 , 39 , 35 , 30 , 24 , 17 , 9 , 0};

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
    case 3: polyms_Ord_3(N,pr,ps,pFmask); break;
    case 4: polyms_Ord_4(N,pr,ps,pFmask); break;
    case 5: polyms_Ord_5(N,pr,ps,pFmask); break;
    case 6: polyms_Ord_6(N,pr,ps,pFmask); break;
    case 7: polyms_Ord_7(N,pr,ps,pFmask); break;
    case 8: polyms_Ord_8(N,pr,ps,pFmask); break;
    default:
        WxLogger *log = WxLogger::get("apollo-root.console");
        WxLogStream infStrm = log->getInfoStream();
        infStrm << "** The polynomial order selected " << N << " is not valid. Max N should be 8." << std::endl;
        break;
    }
}

#endif // WXPNODALDGFUNCTIONS_H
