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

    /* find element-element connections */
    _EtoV = alloc_2d_c<int>(_Klocal,_Vlocal);
    _EToE = alloc_2d_c<int>(_Klocal,_NfE);
    _EToF = alloc_2d_c<int>(_Klocal,_NfE);
    _ETETF = alloc_2d_c<int>(_Klocal,2*_NfE);
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
    free_2d_c(_ETETF, _Klocal,2*_NfE);
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

            int num = cells[0]==K? cells[1]: cells[0];
            _ETETF[K][2*F] = num;
            const PetscInt *ff, *cc;
            DMPlexGetCone(dm, num, &ff);
            for(PetscInt face=0; face<_NfE; face++)
            {
                DMPlexGetSupport(dm, ff[face], &cc);
                if(cc[0]==K || cc[1]==K)
                    _ETETF[K][2*F+1] = face;
            }

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
    PetscInt coordSize, off;

    DMGetCoordinates(dm, &coordinates);
    DMGetCoordinateSection(dm, &coordSection);

    PetscInt eStart, eEnd;
    DMPlexGetHeightStratum(dm, 0, &eStart, &eEnd);

    VecGetArray(coordinates, &coords);
    for(unsigned K=eStart; K<eEnd; K++)
    {
        DMPlexVecGetClosure(dm, coordSection, coordinates, K, &coordSize, &coords);
        //PetscSectionGetOffset(coordSection, K, &off);
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
WxpDGGeometry<REAL>::LIFT_flux(REAL *nflux, REAL *nFrhs, REAL *norms)
{
    for(unsigned K=0; K<_NpE*_meqn; K++)
        nflux[K] = 0.0;

    for(unsigned K=0; K<_NpE; K++)
        for(unsigned f1=0; f1<_NfE; f1++)
            for(unsigned f2=0; f2<_NpF; f2++)
                for(unsigned cp=0; cp<_meqn; cp++)
                    nflux[K*_meqn+cp] += _LIFT[f1*_NpF+f2]*nFrhs[(f1*_NpF+f2)*
                            _meqn+cp]/norms[_NfE*f1+2];
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

    for(unsigned comp=0; comp<_meqn; comp++){
        for(unsigned nk=0; nk<_NpE; nk++){
            for(unsigned mk=0; mk<_NpE; mk++){
                DxnDy[nk*_meqn+comp] += rx[nk]*_Dr[nk*_NpE+mk]*Fflux[mk]+
                                        sx[nk]*_Ds[nk*_NpE+mk]*Fflux[mk]+
                                        ry[nk]*_Dr[nk*_NpE+mk]*Gflux[mk]+
                                        sy[nk]*_Ds[nk*_NpE+mk]*Gflux[mk];
            }
        }
    }
}

// instantiations
template class WxpDGGeometry<float>;
template class WxpDGGeometry<double>;
