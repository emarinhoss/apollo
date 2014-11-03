#include "wxpdggeometry.h"

// WarpX lib includes
#include <wxmath.h>
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

    /* find element-element connections */
    _EtoV = alloc_2d_c<int>(_Klocal,_Vlocal);
    _EToE = alloc_2d_c<int>(_Klocal,_NfE);
    _EToF = alloc_2d_c<int>(_Klocal,_NfE);
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
    free_2d_c(_EToE, _Klocal, _NfE);
    free_2d_c(_EToF, _Klocal, _NfE);
    free_2d_c(_xcoord, _Klocal, _NpE);
    free_2d_c(_ycoord, _Klocal, _NpE);
}

template <typename REAL>
void
WxpDGGeometry<REAL>::FacePair2d(DM dm)
{
    PetscInt eStart, eEnd;
    DMPlexGetHeightStratum(_dm, 0, &eStart, &eEnd);
    for(PetscInt K=eStart; K<eEnd; K++)
    {
        const PetscInt *faces, *cells, *vertex;
        DMPlexGetCone(dm, K, &faces);
        for(PetscInt F=0; F<_NfE; F++)
        {
            _EToF[K][F] = faces[F];
            DMPlexGetSupport(dm, faces[F], &cells);
            _EToE[K][F] = cells[0]==K? cells[1]: cells[0];
            DMPlexGetCone(dm, faces[F], &vertex);
            _EtoV[K][F] = vertex[0];
        }
    }
}

template <typename REAL>
void
WxpDGGeometry<REAL>::CalculateNodeCoordinates2d(DM dm)
{
    Vec coordinates;
    PetscSection coordSection;
    PetscScalar *coords;
    PetscInt coordSize;

    DMGetCoordinates(dm, &coordinates);
    DMGetCoordinateSection(dm, &coordSection);

    PetscInt eStart, eEnd;
    DMPlexGetHeightStratum(dm, 0, &eStart, &eEnd);

    VecGetArray(coordinates, &coords);
    for(unsigned K=eStart; K<eEnd; K++)
    {
        DMPlexVecGetClosure(dm, coordSection, coordinates, K, &coordSize, &coords);
        // coords is returned as coords[x1,y1,x2,y2,x3,y3]
        for(unsigned node=0; node<_NpE; node++)
        {
            REAL r = _r[node];
            REAL s = _s[node];

            _xcoord[K][node] = 0.5*(-coords[0]*(r+s) + coords[2]*(1.+r) + coords[4]*(1.+ s));
            _ycoord[K][node] = 0.5*(-coords[1]*(r+s) + coords[3]*(1.+r) + coords[5]*(1.+ s));
        }
        //DMPlexVecRestoreClosure(dm, coordSection, coordinates, K, &coordSize, &coords);
    }
    VecRestoreArray(coordinates, &coords);
}

template <typename REAL>
void
WxpDGGeometry<REAL>::GeometricFactors2d(int k, REAL *drdx, REAL *dsdx, REAL *drdy, REAL *dsdy, REAL *J)
{
    REAL x1 = _xcoord[k][0], y1 =  _ycoord[k][0];
    REAL x2 = _xcoord[k][1], y2 =  _ycoord[k][1];
    REAL x3 = _xcoord[k][2], y3 =  _ycoord[k][2];

    REAL dxdr = (x2-x1)/2,  dxds = (x3-x1)/2;
    REAL dydr = (y2-y1)/2,  dyds = (y3-y1)/2;

    /* Jacobian of coordinate mapping */
    *J = -dxds*dydr + dxdr*dyds;

    if(*J<0)
      printf("warning: J = %lg\n", *J);

    /* inverted Jacobian matrix for coordinate mapping */
    *drdx =  dyds/(*J);
    *dsdx = -dydr/(*J);
    *drdy = -dxds/(*J);
    *dsdy =  dxdr/(*J);
}

template <typename REAL>
void
WxpDGGeometry<REAL>::Normals2d(int k, REAL *nx, REAL *ny, REAL *sJ)
{
    int f;

    REAL x1 = _xcoord[k][0], y1 = _ycoord[k][0];
    REAL x2 = _xcoord[k][1], y2 = _ycoord[k][1];
    REAL x3 = _xcoord[k][2], y3 = _ycoord[k][2];

    nx[0] =  (y2-y1);  ny[0] = -(x2-x1);
    nx[1] =  (y3-y2);  ny[1] = -(x3-x2);
    nx[2] =  (y1-y3);  ny[2] = -(x1-x3);

    for(f=0;f<_NfE;++f)
    {
      sJ[f] = sqrt(nx[f]*nx[f]+ny[f]*ny[f]);
      nx[f] /= sJ[f];
      ny[f] /= sJ[f];
      sJ[f] /= 2.;
    }
}

// instantiations
template class WxpDGGeometry<float>;
template class WxpDGGeometry<double>;
