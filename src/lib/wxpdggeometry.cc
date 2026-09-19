#include "wxpdggeometry.h"

// WarpX lib includes
#include <wxmath.h>
#include <wxlogger.h>
#include <wxlogstream.h>
#include "wxpnodaldgfunctions.h"
#include "petsc_compat.h"  // PETSc API compatibility for version 3.19+

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

    _Dr   = alloc_1d<REAL>(_NpE*_NpE);
    _Ds   = alloc_1d<REAL>(_NpE*_NpE);
    _Drw  = alloc_1d<REAL>(_NpE*_NpE);
    _Dsw  = alloc_1d<REAL>(_NpE*_NpE);
    _Vand = alloc_1d<REAL>(_NpE*_NpE);
    _IVand= alloc_1d<REAL>(_NpE*_NpE);
    _LIFT = alloc_1d<REAL>(_NpE*_NpF*_NfE);
    _Fmask= alloc_1d<int>(_NfE*_NpF);

    nodalDGfunctions(_SpOr, _r, _s, _Dr, _Ds, _Drw, _Dsw, _Vand, _IVand, _LIFT, _Fmask);
    //
    _rmin = fabs(_r[1]-_r[0]);

    PetscInt eStart, eEnd, vStart, vEnd;
    DMPlexGetHeightStratum(_dm, 0, &eStart, &eEnd);
    DMPlexGetDepthStratum(_dm, 0, &vStart, &vEnd);
    _Klocal = eEnd - eStart;
    _Vlocal = vEnd - vStart;

    /* find element to element connections */
    _EtoV   = alloc_2d_c<int>(_Klocal,_Vlocal);
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
    delete [] _Drw;
    delete [] _Dsw;
    delete [] _Vand;
    delete [] _IVand;
    delete [] _LIFT;
    delete [] _Fmask;
    free_2d_c(_EtoV, _Klocal, _Vlocal);
    free_2d_c(_ETETF, _Klocal,2*_NfE);
    free_2d_c(_xcoord, _Klocal, _NpE);
    free_2d_c(_ycoord, _Klocal, _NpE);
}

template <typename REAL>
void
WxpDGGeometry<REAL>::FacePair2d(DM dm)
{
    Mat FtoV, FtoF;
    DM unint;
    // Build the Element connectivity matrix, EtoV
    PetscInt eStart, eEnd, eEndInt;
    DMPlexUninterpolate(dm,&unint);
    DMPlexGetHeightStratum(dm, 0, &eStart, &eEnd);
    DMPlexGetHybridBounds(dm, &eEndInt, NULL, NULL, NULL);
    for(PetscInt K=eStart; K<eEnd; K++)
    {
        const PetscInt *vertex;
        DMPlexGetCone(unint, K, &vertex);
        for(unsigned vert=0; vert<3; vert++){
            _EtoV[K][vert] = vertex[vert]-_Klocal+1;}
    }

    /** ===================================== */
    /** Build Face to Vertex connection, FtoV */
    /** ===================================== */

    MatCreateSeqAIJ(PETSC_COMM_SELF,_NfE*_Klocal,_Vlocal,2,PETSC_NULL,&FtoV);

    // zero values
    MatZeroEntries(FtoV);

    // Build connection
    int vn[3][2] = {{0,1},{1,2},{2,0}};
    int sk = 0;
    for(unsigned elem=0; elem<_Klocal; elem++)
        for(unsigned face=0; face<_NfE; face++)
        {
            for(unsigned node=0; node<2; node++)
                MatSetValue(FtoV, sk, _EtoV[elem][vn[face][node]]-1, 1, INSERT_VALUES);
            sk++;
        }
    MatAssemblyBegin(FtoV, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(FtoV, MAT_FINAL_ASSEMBLY);

    /** ===================================== */
    /** Build Face to Face connection, FtoF */
    /** ===================================== */

    MatMatTransposeMult(FtoV,FtoV,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&FtoF);

    // substract diagonal contribution
    for(unsigned K1=0; K1<_NfE*_Klocal; K1++)
        MatSetValue(FtoF, K1, K1, -2., ADD_VALUES);
    MatAssemblyBegin(FtoF, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(FtoF, MAT_FINAL_ASSEMBLY);

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
     * If the values of e1 and f1 are negative (-1), this face connects to
     * physical boundary and a boundary condition must be applied at this face.
     * ===================================== */

    PetscInt ncols;
    const PetscInt    *cols;
    const PetscScalar *vals;
    //
    int f1[_NfE*_Klocal], f2[_NfE*_Klocal];
    _totNFace = 0;
    for(unsigned K1=0; K1<_NfE*_Klocal; K1++)
    {
        MatGetRow(FtoF,K1,&ncols,&cols,&vals);
        for(unsigned K2=0; K2<ncols; K2++)
            if(vals[K2]==2.)
            {
                f1[_totNFace] = K1;
                f2[_totNFace++] = cols[K2];
            }
        MatRestoreRow(FtoF,K1,&ncols,&cols,&vals);
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

    MatDestroy(&FtoV);
    MatDestroy(&FtoF);
}

template <typename REAL>
void
WxpDGGeometry<REAL>::CalculateNodeCoordinates2d(DM dm)
{
    Vec coordinates;
    DM unint;
    PetscSection coordSection, defaultSec;
    PetscScalar *coords;
    const PetscInt *pcone;
    //PetscInt coordSize;
    _dtscale = 1.0e6;

    DMGetCoordinatesLocal(dm, &coordinates);
    DMGetCoordinateSection(dm, &coordSection);
    DMGetDefaultSection(dm, &defaultSec);

    PetscInt eStart, eEnd, eEndInt;
    DMPlexGetHeightStratum(dm, 0, &eStart, &eEnd);
    DMPlexGetHybridBounds(dm, &eEndInt, NULL, NULL, NULL);
    DMPlexUninterpolate(dm,&unint);

    VecGetArray(coordinates, &coords);
    for(unsigned K=eStart; K<eEnd; K++)
    {
        //DMPlexVecGetClosure(dm, coordSection, coordinates, K, &coordSize, &coords);
        //PetscSectionGetOffset(defaultSec, K, &off);
        // coords is returned as coords[x1,y1,x2,y2,x3,y3]
        DMPlexGetCone(unint,K,&pcone);
        REAL p1x = coords[2*(pcone[0]-eEnd)];
        REAL p1y = coords[2*(pcone[0]-eEnd)+1];
        REAL p2x = coords[2*(pcone[1]-eEnd)];
        REAL p2y = coords[2*(pcone[1]-eEnd)+1];
        REAL p3x = coords[2*(pcone[2]-eEnd)];
        REAL p3y = coords[2*(pcone[2]-eEnd)+1];

        for(unsigned node=0; node<_NpE; node++)
        {
            REAL r = _r[node];
            REAL s = _s[node];

            _xcoord[K][node] = 0.5*(-p1x*(r+s) + p2x*(1.+r) + p3x*(1.+ s));
            _ycoord[K][node] = 0.5*(-p1y*(r+s) + p2y*(1.+r) + p3y*(1.+ s));
        }
        //DMPlexVecRestoreClosure(dm, coordSection, coordinates, K, &coordSize, &coords);
        REAL len1 = sqrt(pow(p1x-p2x,2)+pow(p1y-p2y,2));
        REAL len2 = sqrt(pow(p2x-p3x,2)+pow(p2y-p3y,2));
        REAL len3 = sqrt(pow(p3x-p1x,2)+pow(p3y-p1y,2));
        REAL sper = 0.5*(len1+len2+len3);
        REAL Area = sqrt(fabs(sper*(sper-len1)*(sper-len2)*(sper-len3)));

        // Compute minimum scale using radius of inscribed circle
        _dtscale = dmin(_dtscale,Area/sper);
    }
    //VecRestoreArray(coordinates, &coords);

}

//template <typename REAL>
//void
//WxpDGGeometry<REAL>::GeometricFactors2d(int k, REAL geom[])
//{
//    REAL x1 = _xcoord[k][_Fmask[0*_NpF]], y1 =  _ycoord[k][_Fmask[0*_NpF]];
//    REAL x2 = _xcoord[k][_Fmask[1*_NpF]], y2 =  _ycoord[k][_Fmask[1*_NpF]];
//    REAL x3 = _xcoord[k][_Fmask[2*_NpF]], y3 =  _ycoord[k][_Fmask[2*_NpF]];

//    REAL dxdr = (x2-x1)/2,  dxds = (x3-x1)/2;
//    REAL dydr = (y2-y1)/2,  dyds = (y3-y1)/2;

//    /* Jacobian of coordinate mapping */
//    REAL J = -dxds*dydr + dxdr*dyds;

//    if(J<=0){
//        WxLogger *l = WxLogger::get("apollo-root.console");
//        WxLogStream errStrm = l->getErrorStream();
//        errStrm << "Error: Jacobian determinant for element " << k << " is " << J;
//        exit(1); // abort execution
//    }

//    /* inverted Jacobian matrix for coordinate mapping */
//    geom[0] =  dyds/(J);
//    geom[1] = -dydr/(J);
//    geom[2] = -dxds/(J);
//    geom[3] =  dxdr/(J);
//    geom[4] =  J;
//}

//template <typename REAL>
//void
//WxpDGGeometry<REAL>::GeomFacs2d(int k, REAL *rx, REAL *sx, REAL *ry, REAL *sy, REAL *J)
//{
//    REAL x[_NpE], y[_NpE];
//    REAL xr[_NpE], yr[_NpE], xs[_NpE], ys[_NpE];

//    for(unsigned n=0; n<_NpE; n++){
//        x[n] = _xcoord[k][n];
//        y[n] = _ycoord[k][n];
//        xr[n]=0.0; yr[n]=0.0; xs[n]=0.0; ys[n]=0.0; J[n]=0.0;
//        rx[n]=0.0; ry[n]=0.0; sx[n]=0.0; sy[n]=0.0;
//    }

//    // Do matrix vector products
//    for(unsigned m=0; m<_NpE; m++)
//        for(unsigned n=0; n<_NpE; n++){
//            xr[m] += _Dr[m*_NpE+n]*x[n];
//            xs[m] += _Ds[m*_NpE+n]*x[n];
//            yr[m] += _Dr[m*_NpE+n]*y[n];
//            ys[m] += _Ds[m*_NpE+n]*y[n];
//        }

//    for(unsigned m=0; m<_NpE; m++){
//        J[m] = -xs[m]*yr[m]+xr[m]*ys[m];
//        rx[m]=  ys[m]/J[m];
//        sx[m]= -yr[m]/J[m];
//        ry[m]= -xs[m]/J[m];
//        sy[m]=  xr[m]/J[m];
//    }
//}

//template <typename REAL>
//void
//WxpDGGeometry<REAL>::Normals2d(int k, REAL norms[])
//{
//    int f;

//    REAL x1 = _xcoord[k][_Fmask[0*_NpF]], y1 = _ycoord[k][_Fmask[0*_NpF]];
//    REAL x2 = _xcoord[k][_Fmask[1*_NpF]], y2 = _ycoord[k][_Fmask[1*_NpF]];
//    REAL x3 = _xcoord[k][_Fmask[2*_NpF]], y3 = _ycoord[k][_Fmask[2*_NpF]];

//    norms[0] =  (y2-y1);  norms[1] = -(x2-x1);
//    norms[3] =  (y3-y2);  norms[4] = -(x3-x2);
//    norms[6] =  (y1-y3);  norms[7] = -(x1-x3);

//    for(f=0;f<_NfE;++f)
//    {
//      REAL sJ = sqrt(norms[_NfE*f]*norms[_NfE*f]+norms[_NfE*f+1]*norms[_NfE*f+1]);
//      if(sJ<=0){
//          WxLogger *l = WxLogger::get("apollo-root.console");
//          WxLogStream errStrm = l->getErrorStream();
//          errStrm << "Error: Edge length of " << sJ << " found in element " << k;
//          exit(1); // abort execution
//      }
//      norms[_NfE*f]   /= sJ;
//      norms[_NfE*f+1] /= sJ;
//      norms[_NfE*f+2] = sJ/2.;
//    }
//}

template <typename REAL>
void
WxpDGGeometry<REAL>::FaceNodesNormals2d(int k, REAL *nx, REAL *ny, REAL *Fscale)
{
    REAL x[_NpE], y[_NpE];
    REAL xr[_NpE], yr[_NpE], xs[_NpE], ys[_NpE], J[_NpE], sJ[_NpF*_NfE];
    int Fmask[_NpF*_NfE];

    returnFmask(Fmask);

    for(unsigned n=0; n<_NpE; n++){
        x[n] = _xcoord[k][n];
        y[n] = _ycoord[k][n];
        xr[n]=0.0; yr[n]=0.0; xs[n]=0.0; ys[n]=0.0; J[n]=0.0;
    }

    // Do matrix vector products
    for(unsigned m=0; m<_NpE; m++)
        for(unsigned n=0; n<_NpE; n++){
            xr[m] += _Dr[m*_NpE+n]*x[n];
            xs[m] += _Ds[m*_NpE+n]*x[n];
            yr[m] += _Dr[m*_NpE+n]*y[n];
            ys[m] += _Ds[m*_NpE+n]*y[n];
        }

    for(unsigned m=0; m<_NpE; m++)
    {
            J[m] = -xs[m]*yr[m]+xr[m]*ys[m];
            if(J[m]<=0){
                WxLogger *l = WxLogger::get("apollo-root.console");
                WxLogStream errStrm = l->getErrorStream();
                errStrm << "Error: Jacobian determinant for element " << k
                        << " node " << m << " is " << J[m]
                        << ". A zero Jacobian across every element usually means "
                           "the differentiation operators are zero rather than "
                           "that the mesh is degenerate; a negative one means an "
                           "inverted element in the mesh.";
                exit(1); // abort execution
            }
    }

    // Face 1
    for(unsigned fid=0; fid<_NpF; fid++)
    {
        nx[fid] = yr[Fmask[fid]];
        ny[fid] = -xr[Fmask[fid]];
    }

    // Face 2
    for(unsigned fid=_NpF; fid<2*_NpF; fid++)
    {
        nx[fid] = ys[Fmask[fid]]-yr[Fmask[fid]];
        ny[fid] = -xs[Fmask[fid]]+xr[Fmask[fid]];
    }

    // Face 3
    for(unsigned fid=2*_NpF; fid<3*_NpF; fid++)
    {
        nx[fid] = -ys[Fmask[fid]];
        ny[fid] = xs[Fmask[fid]];
    }

    for(unsigned fid=0; fid<_NpF*_NfE; fid++)
    {
        sJ[fid]  = sqrt(nx[fid]*nx[fid]+ny[fid]*ny[fid]);
        if(sJ[fid]<=0){
                  WxLogger *l = WxLogger::get("apollo-root.console");
                  WxLogStream errStrm = l->getErrorStream();
                  errStrm << "Error: Edge length of " << sJ[fid] << " found in element " << k;
                  exit(1); // abort execution
              }
        nx[fid] /= sJ[fid];
        ny[fid] /= sJ[fid];
        Fscale[fid] = sJ[fid]/J[Fmask[fid]];
    }

}

template <typename REAL>
void
WxpDGGeometry<REAL>::LIFT_flux(int k, REAL *nflux, REAL *nFrhs)
{
    for(unsigned K=0; K<_NpE*_meqn; K++)
        nflux[K] = 0.0;

    for(unsigned nodes=0; nodes<_NpE; nodes++)
        for(unsigned faceNode=0; faceNode<_NpF*_NfE; faceNode++)
            for(unsigned comp=0; comp<_meqn; comp++)
                nflux[nodes*_meqn+comp] += _LIFT[nodes*_NpF*_NfE+faceNode]*nFrhs[faceNode*_meqn+comp];
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
        rx[n]=0.0; ry[n]=0.0; sx[n]=0.0; sy[n]=0.0;}

    for(unsigned nk=0; nk<_NpE*_meqn; nk++)
        DxnDy[nk] = 0.0;

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
                DxnDy[nk*_meqn+comp] += rx[nk]*_Dr[nk*_NpE+mk]*Fflux[mk*_meqn+comp]+
                                        sx[nk]*_Ds[nk*_NpE+mk]*Fflux[mk*_meqn+comp]+
                                        ry[nk]*_Dr[nk*_NpE+mk]*Gflux[mk*_meqn+comp]+
                                        sy[nk]*_Ds[nk*_NpE+mk]*Gflux[mk*_meqn+comp];

}

template <typename REAL>
void
WxpDGGeometry<REAL>::calculateFilter(REAL *filter, REAL *filterMatrix)
{
    REAL fitTimesIVand[_NpE*_NpE], filterIn[_NpE*_NpE];

    // zero entries
    for(unsigned k=0; k<_NpE*_NpE; k++){
        fitTimesIVand[k] = 0.0;
        filterMatrix[k] = 0.0;
        filterIn[k] = 0.0;}

    // create a filter diagonal matrix
    for(unsigned k=0; k<_NpE; k++)
        filterIn[k*_NpE+k] = filter[k];

    // do Vand*diag(filter)*IVand
    for(unsigned row=0; row<_NpE; row++)
        for(unsigned col=0;col<_NpE;col++)
            for(unsigned k=0; k<_NpE; k++)
                fitTimesIVand[row*_NpE+col] += filterIn[row*_NpE+k]*_IVand[col+k*_NpE];

    for(unsigned col=0; col<_NpE; col++)
        for(unsigned row=0;row<_NpE;row++)
            for(unsigned k=0; k<_NpE; k++)
                filterMatrix[row*_NpE+col] += _Vand[row*_NpE+k]*fitTimesIVand[col+k*_NpE];
}

// instantiations
//template class WxpDGGeometry<float>;
template class WxpDGGeometry<double>;
