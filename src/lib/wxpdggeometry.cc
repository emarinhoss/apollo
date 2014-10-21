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

    p_Dr = alloc_2d_c<REAL>(_NpE,_NpE);
    p_Ds = alloc_2d_c<REAL>(_NpE,_NpE);
    p_LIFT = alloc_2d_c<REAL>(_NpE,_NpF*_NfE);
    p_Fmask = alloc_2d_c<int>(_NfE,_NpF);

    //matrices
    _Dr = matrix_REAL(_NpE,_NpE);
    _Ds = matrix_REAL(_NpE,_NpE);
    _LIFT = matrix_REAL(_NpE,_NpF*_NfE);
    _Fmask = matrix_INT(_NfE,_NpF);

    nodalDGfunctions(_SpOr, _r, _s, p_Dr, p_Ds, p_LIFT, p_Fmask);

    for(unsigned i=0; i<_NpE; i++){
        for(unsigned j=0; j<_NpE; j++){
            _Dr(i,j) = p_Dr[i][j];
            _Ds(i,j) = p_Ds[i][j];}

        for(unsigned k=0; k<_NpF*_NfE; k++)
            _LIFT(i,k) = p_LIFT[i][k];
    }

    for(unsigned i=0; i<_NfE; i++)
        for(unsigned j=0; j<_NpF; j++)
            _Fmask(i,j) = p_Fmask[i][j];

    PetscInt eStart, eEnd, vStart, vEnd;
    DMPlexGetHeightStratum(_dm, 0, &eStart, &eEnd);
    DMPlexGetDepthStratum(_dm, 0, &vStart, &vEnd);
    _Klocal = eEnd - eStart;
    _Vlocal = vStart - vEnd;

    /* find element-element connections */
    _EtoV = alloc_2d_c<unsigned>(_Klocal,_Vlocal);
    _EToE = alloc_2d_c<unsigned>(_Klocal,_NfE);
    _EToF = alloc_2d_c<unsigned>(_Klocal,_NfE);
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
    free_2d_c(_EtoV, _Klocal, _Vlocal);
    free_2d_c(_EToE, _Klocal, _NfE);
    free_2d_c(_EToF, _Klocal, _NfE);
    free_2d_c(_xcoord, _Klocal, _NpE);
    free_2d_c(_ycoord, _Klocal, _NpE);
    free_2d_c(p_Dr, _NpE, _NpE);
    free_2d_c(p_Ds, _NpE, _NpE);
    free_2d_c(p_LIFT, _NpE, _NpF*_NfE);
    free_2d_c(p_Fmask, _NfE, _NpF);
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
    PetscScalar *coords;

    DMGetCoordinatesLocal(dm, &coordinates);
    int _dim = 2;

    PetscInt eStart, eEnd;
    DMPlexGetHeightStratum(dm, 0, &eStart, &eEnd);

    VecGetArray(coordinates, &coords);
    for(unsigned K=eStart; K<eEnd; K++)
    {
        for(unsigned node=0; node<_NpE; node++)
        {
            REAL r = _r[node];
            REAL s = _s[node];
            REAL Gx1 = coords[_EtoV[K][0]*_dim];
            REAL Gx2 = coords[_EtoV[K][1]*_dim];
            REAL Gx3 = coords[_EtoV[K][2]*_dim];
            REAL Gy1 = coords[_EtoV[K][0]*_dim+1];
            REAL Gy2 = coords[_EtoV[K][1]*_dim+1];
            REAL Gy3 = coords[_EtoV[K][2]*_dim+1];

            _xcoord[K][node] = 0.5*(-Gx1*(r+s) + Gx2*(1.+r) + Gx3*(1.+ s));
            _ycoord[K][node] = 0.5*(-Gy1*(r+s) + Gy2*(1.+r) + Gy3*(1.+ s));
        }
    }
    VecRestoreArray(coordinates, &coords);
}

// instantiations
template class WxpDGGeometry<float>;
template class WxpDGGeometry<double>;
