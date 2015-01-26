#include "wxpdggeometry.h"

// WarpX lib includes
#include <wxmath.h>
#include <wxlogger.h>
#include <wxlogstream.h>
#include "wxpnodaldgfunctions.h"

template <typename REAL>
WxpDGGeometry<REAL>::WxpDGGeometry(DM dm, unsigned meqn, unsigned Spor)
    : _meqn(meqn), _SpOr(Spor), _dm(dm)
{
    _NpE = (_SpOr+1)*(_SpOr+2)/2;
    _NpF = (_SpOr+1);
    _NfE = 3; // for now only triangular elements (will generalize later!)    

    // allocate memory
    _r = alloc_1d<REAL>(_NpE);
    _s = alloc_1d<REAL>(_NpE);

    _Dr = alloc_1d<REAL>(_NpE*_NpE);
    _Ds = alloc_1d<REAL>(_NpE*_NpE);
    _LIFT = alloc_1d<REAL>(_NpE*_NpF*_NfE);
    _Fmask = alloc_1d<int>(_NfE*_NpF);

    nodalDGfunctions(_SpOr, _r, _s, _Dr, _Ds, _LIFT, _Fmask);

    PetscInt eStart, eEnd, vStart, vEnd;
    DMPlexGetHeightStratum(_dm, 0, &eStart, &eEnd);
    DMPlexGetDepthStratum(_dm, 0, &vStart, &vEnd);
    _Klocal = eEnd - eStart;
    _Vlocal = vEnd - vStart;

    /* find element to element connections */
    _EtoV   = alloc_2d_c<int>(_Klocal,_Vlocal);
    _FToV   = alloc_2d_c<int>(_NfE*_Klocal,_Vlocal);
    _FToV_t = alloc_2d_c<int>(_Vlocal,_NfE*_Klocal);
    _FToF   = alloc_2d_c<int>(_NfE*_Klocal,_NfE*_Klocal);
    _ETETF  = alloc_2d_c<int>(_Klocal,2*_NfE);
    FacePair2d(_dm);

    // Find node coordinates
    _xcoord = alloc_2d_c<REAL>(_Klocal,_NpE);
    _ycoord = alloc_2d_c<REAL>(_Klocal,_NpE);
    CalculateNodeCoordinates2d(_dm);

}

template <typename REAL>
WxpDGGeometry<REAL>::~WxpDGGeometry()
{
    delete [] _r;
    delete [] _s;
    delete [] _Dr;
    delete [] _Ds;
    delete [] _LIFT;
    delete [] _Fmask;
    free_2d_c(_EtoV, _Klocal, _Vlocal);
    free_2d_c(_FToV, _NfE*_Klocal, _Vlocal);
    free_2d_c(_FToV_t, _Vlocal, _NfE*_Klocal);
    free_2d_c(_FToF, _NfE*_Klocal, _NfE*_Klocal);
    free_2d_c(_ETETF, _Klocal,2*_NfE);
    free_2d_c(_xcoord, _Klocal, _NpE);
    free_2d_c(_ycoord, _Klocal, _NpE);
}

template <typename REAL>
void
WxpDGGeometry<REAL>::FacePair2d(DM dm)
{
    // Build the Element connectivity matrix, EtoV
    PetscInt eStart, eEnd, eEndInt;
    DMPlexGetHeightStratum(_dm, 0, &eStart, &eEnd);
    DMPlexGetHybridBounds(dm, &eEndInt, NULL, NULL, NULL);
    for(PetscInt K=eStart; K<eEndInt; K++)
    {
        const PetscInt *vertex;
        DMPlexGetCone(dm, K, &vertex);
        for(unsigned vert=0; vert<3; vert++)
            _EtoV[K][vert] = vertex[vert]-_Klocal+1;
    }

    /** ===================================== */
    /** Build Face to Vertex connection, FtoV */
    /** ===================================== */

    // zero values
    for(int K=0; K<_Klocal*_NfE; K++)
        for(int V=0; V<_Vlocal; V++)
        {
            _FToV[K][V] = 0;
            _FToV_t[V][K] = 0;
        }

    // Build connection
    int vn[3][2] = {{0,1},{1,2},{2,0}};
    int sk = 0;
    for(unsigned elem=0; elem<_Klocal; elem++)
        for(unsigned face=0; face<_NfE; face++)
        {
            for(unsigned node=0; node<2; node++)
            {
                _FToV[sk][_EtoV[elem][vn[face][node]]-1]   = 1;
                _FToV_t[_EtoV[elem][vn[face][node]]-1][sk] = 1;
            }
            sk++;
        }

    /** ===================================== */
    /** Build Face to Face connection, FtoF */
    /** ===================================== */

    // zero values
    for(unsigned K=0; K<_NfE*_Klocal; K++)
        for(unsigned V=0; V<_NfE*_Klocal; V++)
            _FToF[K][V] = 0;

    // Build connection
    for(int K1=0; K1<_NfE*_Klocal; K1++)
        for(int K2=0; K2<_NfE*_Klocal; K2++)
            for(int V1=0; V1<_Vlocal; V1++)
                _FToF[K1][K2] += _FToV[K1][V1]*_FToV_t[V1][K2];

    // substract diagonal contribution
    for(unsigned K1=0; K1<_NfE*_Klocal; K1++)
        _FToF[K1][K1] += -2;

    /** =====================================
     * Build Element to Element to Face connection, _ETETF.
     * The rows represent a specific element K, the two first columns
     * tell information about the first face (zeroth face), the next 2
     * about the next face and so on ...
     * For the zeroth face, the first number tells the element on the outside
     * of the face, and the second number represents the face number of the outside
     * element.
     * Element |  Face0  |  Face1  |  Face2  |
     * ---------------------------------
     *     k   | e1 | f1 | ....
     *
     * Face 0 of element k is connect to element e1 at face f1 (of element e1).
     * If the values of values of e1 and f1 are negative (-1), this face connects to
     * physical boundary and a boundary condition must be applied at this face.
     * ===================================== */

    //
    int f1[_NfE*_Klocal], f2[_NfE*_Klocal];
    _totNFace = 0;
    for(unsigned K1=0; K1<_NfE*_Klocal; K1++)
        for(unsigned K2=0; K2<_NfE*_Klocal; K2++)
            if(_FToF[K1][K2]==2)
            {
                f1[_totNFace] = K1;
                f2[_totNFace++] = K2;
            }

    int elem1[_totNFace], elem2[_totNFace], face1[_totNFace], face2[_totNFace];
    for(unsigned face=0; face<_totNFace; face++)
    {
        elem1[face] = floor(f1[face]/_NfE);
        elem2[face] = floor(f2[face]/_NfE);

        face1[face] = f1[face]%_NfE;
        face2[face] = f2[face]%_NfE;
    }

    // Make all values -1. Only the faces not at physical boundaries are changed.
    for(unsigned k1=0; k1<_Klocal; k1++)
        for(unsigned f1=0; f1<2*_NfE; f1++)
        {
            _ETETF[k1][f1] = -1;
        }

    // assign values
    for(unsigned kk=0; kk<_totNFace; kk++)
    {
        _ETETF[elem1[kk]][2*face1[kk]]   = elem2[kk];
        _ETETF[elem1[kk]][2*face1[kk]+1] = face2[kk];
    }

}

template <typename REAL>
void
WxpDGGeometry<REAL>::CalculateNodeCoordinates2d(DM dm)
{
    Vec coordinates;
    PetscSection coordSection, defaultSec;
    PetscScalar *coords;
    const PetscInt *pcone;
    //PetscInt coordSize;

    DMGetCoordinatesLocal(dm, &coordinates);
    DMGetCoordinateSection(dm, &coordSection);
    DMGetDefaultSection(dm, &defaultSec);

    PetscInt eStart, eEnd, eEndInt;
    DMPlexGetHeightStratum(dm, 0, &eStart, &eEnd);
    DMPlexGetHybridBounds(dm, &eEndInt, NULL, NULL, NULL);

    VecGetArray(coordinates, &coords);
    for(unsigned K=eStart; K<eEndInt; K++)
    {
        //DMPlexVecGetClosure(dm, coordSection, coordinates, K, &coordSize, &coords);
        //PetscSectionGetOffset(defaultSec, K, &off);
        // coords is returned as coords[x1,y1,x2,y2,x3,y3]
        DMPlexGetCone(dm,K,&pcone);
        REAL p1x = coords[2*(pcone[0]-eEndInt)];
        REAL p1y = coords[2*(pcone[0]-eEndInt)+1];
        REAL p2x = coords[2*(pcone[1]-eEndInt)];
        REAL p2y = coords[2*(pcone[1]-eEndInt)+1];
        REAL p3x = coords[2*(pcone[2]-eEndInt)];
        REAL p3y = coords[2*(pcone[2]-eEndInt)+1];

        for(unsigned node=0; node<_NpE; node++)
        {
            REAL r = _r[node];
            REAL s = _s[node];

            _xcoord[K][node] = 0.5*(-p1x*(r+s) + p2x*(1.+r) + p3x*(1.+ s));
            _ycoord[K][node] = 0.5*(-p1y*(r+s) + p2y*(1.+r) + p3y*(1.+ s));
        }
        //DMPlexVecRestoreClosure(dm, coordSection, coordinates, K, &coordSize, &coords);
    }
    //VecRestoreArray(coordinates, &coords);
}

template <typename REAL>
void
WxpDGGeometry<REAL>::GeometricFactors2d(int k, REAL geom[])
{
    REAL x1 = _xcoord[k][0], y1 =  _ycoord[k][0];
    REAL x2 = _xcoord[k][1], y2 =  _ycoord[k][1];
    REAL x3 = _xcoord[k][2], y3 =  _ycoord[k][2];

    REAL dxdr = (x2-x1)/2,  dxds = (x3-x1)/2;
    REAL dydr = (y2-y1)/2,  dyds = (y3-y1)/2;

    /* Jacobian of coordinate mapping */
    REAL J = -dxds*dydr + dxdr*dyds;

    if(J<=0){
        WxLogger *l = WxLogger::get("apollo-root.console");
        WxLogStream errStrm = l->getErrorStream();
        errStrm << "Error: Jacobian determinant for element " << k << " is " << J;
        exit(1); // abort execution
    }

    /* inverted Jacobian matrix for coordinate mapping */
    geom[0] =  dyds/(J);
    geom[1] = -dydr/(J);
    geom[2] = -dxds/(J);
    geom[3] =  dxdr/(J);
    geom[4] =  J;
}

template <typename REAL>
void
WxpDGGeometry<REAL>::Normals2d(int k, REAL norms[])
{
    int f;

    REAL x1 = _xcoord[k][0], y1 = _ycoord[k][0];
    REAL x2 = _xcoord[k][1], y2 = _ycoord[k][1];
    REAL x3 = _xcoord[k][2], y3 = _ycoord[k][2];

    norms[0] =  (y2-y1);  norms[1] = -(x2-x1);
    norms[3] =  (y3-y2);  norms[4] = -(x3-x2);
    norms[6] =  (y1-y3);  norms[7] = -(x1-x3);

    for(f=0;f<_NfE;++f)
    {
      REAL sJ = sqrt(norms[_NfE*f]*norms[_NfE*f]+norms[_NfE*f+1]*norms[_NfE*f+1]);
      if(sJ<=0){
          WxLogger *l = WxLogger::get("apollo-root.console");
          WxLogStream errStrm = l->getErrorStream();
          errStrm << "Error: Edge length of " << sJ << " found in element " << k;
          exit(1); // abort execution
      }
      norms[_NfE*f]   /= sJ;
      norms[_NfE*f+1] /= sJ;
      norms[_NfE*f+2] = sJ/2.;
    }
}

template <typename REAL>
void
WxpDGGeometry<REAL>::LIFT_flux(REAL *nflux, REAL *nFrhs, REAL *Fscale)
{
    for(unsigned K=0; K<_NpE*_meqn; K++)
        nflux[K] = 0.0;

//    for(unsigned points=0; points<_NpF; points++)
//        for(unsigned faces=0; faces<_NfE; faces++)
//            for(unsigned comp=0; comp<_meqn; comp++)
//                nFrhs[faces*_NpF*_meqn+faces*_meqn+comp] = nFrhs[faces*_NpF*_meqn+faces*_meqn+comp]/norms[_NfE*faces+2];

//    for(unsigned nodes=0; nodes<_NpE; nodes++)
//        for(unsigned faceNode=0; faceNode<_NpF*_NfE; faceNode++)
//            for(unsigned comp=0; comp<_meqn; comp++)
//                nflux[nodes*_meqn+comp] += _LIFT[nodes*_NpF*_NfE+faceNode]*nFrhs[faceNode*_meqn+comp];

    for(unsigned Epoints=0; Epoints<_NpE; Epoints++) // loop over all element nodes
        for(unsigned faces=0; faces<_NfE; faces++) // loop over faces/edges
            for(unsigned Fpoints=0; Fpoints<_NpF; Fpoints++) // loop over nodes on faces/edges
                for(unsigned cp=0; cp<_meqn; cp++) // loop over components
                    nflux[Epoints*_meqn+cp] += _LIFT[faces*_NpF+Fpoints]*nFrhs[(faces*_NpF+Fpoints)*
                            _meqn+cp]/Fscale[faces];
}

template <typename REAL>
void
WxpDGGeometry<REAL>::weakDericatives(unsigned K, REAL *DxnDy, REAL *Fflux, REAL *Gflux)
{
    REAL x[_NpE], y[_NpE];
    REAL xr[_NpE], yr[_NpE], xs[_NpE], ys[_NpE], J[_NpE];
    REAL rx[_NpE], ry[_NpE], sx[_NpE], sy[_NpE];

    for(unsigned n=0; n<_NpE; n++){
        x[n] = _xcoord[K][n];
        y[n] = _ycoord[K][n];
        xr[n]=0.0; yr[n]=0.0; xs[n]=0.0; ys[n]=0.0; J[n]=0.0;
        rx[n]=0.0; ry[n]=0.0; sx[n]=0.0; sy[n]=0.0;
        DxnDy[n] = 0.0;
    }

    // Do matrix vector products
    for(unsigned m=0; m<_NpE; m++)
        for(unsigned n=0; n<_NpE; n++){
            xr[m] += _Dr[m*_NpE+n]*x[n];
            xs[m] += _Ds[m*_NpE+n]*x[n];
            yr[m] += _Dr[m*_NpE+n]*y[n];
            ys[m] += _Ds[m*_NpE+n]*y[n];
        }

    for(unsigned m=0; m<_NpE; m++){
        J[m] = -xs[m]*yr[m]+xr[m]*ys[m];
        rx[m]=  ys[m]/J[m];
        sx[m]= -yr[m]/J[m];
        ry[m]= -xs[m]/J[m];
        sy[m]=  xr[m]/J[m];
    }

    for(unsigned nk=0; nk<_NpE; nk++)
        for(unsigned mk=0; mk<_NpE; mk++)
            for(unsigned comp=0; comp<_meqn; comp++)
                DxnDy[nk*_meqn+comp] += rx[mk]*_Dr[nk*_NpE+mk]*Fflux[mk*_meqn+comp]+
                                        sx[mk]*_Ds[nk*_NpE+mk]*Fflux[mk*_meqn+comp]+
                                        ry[mk]*_Dr[nk*_NpE+mk]*Gflux[mk*_meqn+comp]+
                                        sy[mk]*_Ds[nk*_NpE+mk]*Gflux[mk*_meqn+comp];

}

// instantiations
template class WxpDGGeometry<float>;
template class WxpDGGeometry<double>;
