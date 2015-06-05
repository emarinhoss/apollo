#include "wxnodaldggeometry2d.h"

// WarpX lib includes
#include <wxmath.h>
#include <wxlogger.h>
#include <wxlogstream.h>
#include "wxpnodaldgfunctions.h"
#include "wxNodalDGMatrices.h"

template <typename REAL>
wxNodalDGgeometry2D<REAL>::wxNodalDGgeometry2D(DM dm, unsigned meqn, unsigned Spor)
    : _meqn(meqn), _polyOr(Spor), _dm(dm)
{

    WxLogger *log = WxLogger::get("apollo-root.console");
    WxLogStream infStrm = log->getInfoStream();

    _NpE = (_polyOr+1)*(_polyOr+2)/2;
    infStrm << "** There are " << _NpE << " nodes per Element **" << std::endl;
    _NpF = (_polyOr+1);
    _NfE = 3; // for now only triangular elements (will generalize later!)

    // allocate memory
    _r = alloc_1d<REAL>(_NpE);
    _s = alloc_1d<REAL>(_NpE);
    _Fmask = alloc_1d<int>(_NfE*(_polyOr+1));
    _Dr   = alloc_1d<REAL>(_NpE*_NpE);
    _Ds   = alloc_1d<REAL>(_NpE*_NpE);
    _Vand = alloc_1d<REAL>(_NpE*_NpE);
    _VVT = alloc_1d<REAL>(_NpE*_NpE);
    _Mass = alloc_1d<REAL>(_NpE*_NpE);

    nodalNaturalCoordinates(_polyOr,_r,_s,_Fmask);

    Mat Dr, Ds, Vand, VVT, Mass;
    MatCreateSeqDense(PETSC_COMM_SELF,_NpE,_NpE,PETSC_NULL,&Dr);
    MatCreateSeqDense(PETSC_COMM_SELF,_NpE,_NpE,PETSC_NULL,&Ds);
    MatCreateSeqDense(PETSC_COMM_SELF,_NpE,_NpE,PETSC_NULL,&Vand);
    MatCreateSeqDense(PETSC_COMM_SELF,_NpE,_NpE,PETSC_NULL,&_IVand);
    MatCreateSeqDense(PETSC_COMM_SELF,_NpE,_NpE,PETSC_NULL,&VVT);
    MatCreateSeqDense(PETSC_COMM_SELF,_NpE,_NpE,PETSC_NULL,&Mass);

    // Populate matrices
    Vandermonde2D(_polyOr,_NpE, _r, _s, &Vand);
    infStrm << "** done -- Creating Vandemonde Matrix. **" << std::endl;
    //MatView(_Vand,PETSC_VIEWER_STDOUT_WORLD);
    this->invertMatrix(Vand,&_IVand);
    infStrm << "** done -- Creating Inverse Vandemonde Matrix. **" << std::endl;
    //MatView(_IVand,PETSC_VIEWER_STDOUT_WORLD);
    DifferentiationMatrices2D(_polyOr,_NpE,_r,_s,_IVand,&Dr,&Ds);
    infStrm << "** done -- Creating Differentiation Matrices. **" << std::endl;
    //MatView(_Dr,PETSC_VIEWER_STDOUT_WORLD);
    //MatView(_Dr,PETSC_VIEWER_STDOUT_WORLD);

    // Calculate inverse of mass matrix
    Mat dummy;
    MatTranspose(Vand,MAT_INITIAL_MATRIX,&dummy);
    MatMatMult(Vand,dummy,MAT_REUSE_MATRIX,PETSC_DEFAULT,&VVT);
    this->invertMatrix(VVT,&Mass);
    infStrm << "** done -- Creating Inverse Mass Matrix. **" << std::endl;
//    MatView(VVT,PETSC_VIEWER_STDOUT_WORLD);

    // Find node coordinates for each element
    PetscInt eStart, eEnd, vStart, vEnd;
    DMPlexGetHeightStratum(_dm, 0, &eStart, &eEnd);
    DMPlexGetDepthStratum(_dm, 0, &vStart, &vEnd);
    _Klocal = eEnd - eStart;
    _Vlocal = vEnd - vStart;
    _xcoord = alloc_2d_c<REAL>(_Klocal,_NpE);
    _ycoord = alloc_2d_c<REAL>(_Klocal,_NpE);
    CalculateNodeCoordinates2d(_dm);
    infStrm << "** done -- Calculating Node Coordinates. **" << std::endl;

    /* find element to element connections */
    _EtoV   = alloc_2d_c<int>(_Klocal,3);
    _ETETF  = alloc_2d_c<int>(_Klocal,2*_NfE);
    FacePair2d(_dm);
    infStrm << "** done -- Creating face-to-face connections. **" << std::endl;

    // minimum distance between two nodes, in natural coordinates
    _rmin = fabs(_r[1]-_r[0]);

    // transfer all the Matrices data into row-major C-arrays
    petscMatTOArray(Ds,_Ds);
    petscMatTOArray(Dr,_Dr);
    petscMatTOArray(Vand,_Vand);
    petscMatTOArray(VVT,_VVT);
    petscMatTOArray(Mass,_Mass);

    MatDestroy(&Ds);
    MatDestroy(&Dr);
    MatDestroy(&Vand);
    MatDestroy(&VVT);
    MatDestroy(&dummy);
    MatDestroy(&Mass);

}

template <typename REAL>
void
wxNodalDGgeometry2D<REAL>::petscMatTOArray(Mat A, REAL *array)
{
    PetscInt mcols, mrows;
    MatGetSize(A, &mrows, &mcols);

    PetscInt ncols;
    const PetscInt    *cols;
    const PetscScalar *vals;

    int sk = 0;
    for(unsigned kk=0; kk<mrows; kk++)
    {
        MatGetRow(A,kk,&ncols,&cols,&vals);
        for(int kx=0; kx<ncols; kx++)
            array[sk++] = vals[kx];
    }
}

template <typename REAL>
wxNodalDGgeometry2D<REAL>::~wxNodalDGgeometry2D()
{
    delete [] _r;
    delete [] _s;
    delete [] _Fmask;
    delete [] _Dr;
    delete [] _Ds;
    delete [] _Vand;
    delete [] _VVT;
    delete [] _Mass;
    free_2d_c(_xcoord, _Klocal, _NpE);
    free_2d_c(_ycoord, _Klocal, _NpE);
    free_2d_c(_EtoV, _Klocal, 3);
    free_2d_c(_ETETF, _Klocal, 2*_NfE);
    MatDestroy(&_IVand);
}

template <typename REAL>
void
wxNodalDGgeometry2D<REAL>::invertMatrix(Mat A, Mat *invA)
{
    Mat inpA, B;
    IS is;
    MatFactorInfo iluinfo;
    PetscInt ncols;
    const PetscInt    *cols;
    const PetscScalar *vals;

    MatDuplicate(A,MAT_COPY_VALUES,&inpA);

    // begin by creating a dense matrix B and fill it with the identity matrix
    MatGetRow(A,0,&ncols,&cols,&vals);
    MatCreateSeqDense(PETSC_COMM_SELF,ncols,ncols,PETSC_NULL,&B);
    for (int k=0; k<ncols;k++)
        MatSetValue(B,k,k,1.0,INSERT_VALUES);
    MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);

    MatGetFactor(A,"petsc",MAT_FACTOR_LU,&inpA);
    MatLUFactorSymbolic(inpA,A,is,is,&iluinfo);
    MatLUFactorNumeric(inpA,A,&iluinfo);
//    MatLUFactor(invA,is,is,&iluinfo);
    // Calculate inverse
    MatMatSolve(inpA,B,*invA);

    MatDestroy(&inpA);
    MatDestroy(&B);
}

template <typename REAL>
void
wxNodalDGgeometry2D<REAL>::CalculateNodeCoordinates2d(DM dm)
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

    PetscInt eStart, eEnd, eEndInterior;
    DMPlexGetHeightStratum(dm, 0, &eStart, &eEnd);
    DMPlexGetHybridBounds(dm, &eEndInterior, NULL, NULL, NULL);
    DMPlexUninterpolate(dm,&unint);

    VecGetArray(coordinates, &coords);
    for(unsigned K=eStart; K<eEndInterior; K++)
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
    VecRestoreArray(coordinates, &coords);
    DMDestroy(&unint);
//    VecDestroy(&coordinates);

}

template <typename REAL>
void
wxNodalDGgeometry2D<REAL>::FacePair2d(DM dm)
{
    WxLogger *l = WxLogger::get("apollo-root.console");
    WxLogStream errStrm = l->getErrorStream();

    Mat FtoV, FtoF;
    DM unint;
    // Build the Element connectivity matrix, EtoV
    PetscInt eStart, eEnd, eEndInterior;
    DMPlexUninterpolate(dm,&unint);
    DMPlexGetHeightStratum(dm, 0, &eStart, &eEnd);
    DMPlexGetHybridBounds(dm, &eEndInterior, NULL, NULL, NULL);
    for(PetscInt K=eStart; K<eEndInterior; K++)
    {
        const PetscInt *vertex;
        DMPlexGetCone(unint, K, &vertex);
        for(unsigned vert=0; vert<3; vert++){
            _EtoV[K][vert] = vertex[vert]-_Klocal+1;}
    }

    DMDestroy(&unint);

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

    PetscInt value;
    PetscInt vStart, vEnd;
    DMPlexGetHeightStratum(_dm, 1, &vStart, &vEnd);
    // Get the cells that support this face
    const PetscInt *cells;

    for(unsigned face=vStart; face<vEnd; face++)
    {
        DMPlexGetLabelValue(dm, "Face Sets", face, &value);
        if(value!=-1){
            DMPlexGetSupport(dm, face, &cells);
            for(unsigned f1=0; f1<2*_NfE; f1++)
                _ETETF[cells[0]][f1] = -value;
        }
//        int AAA = 0;
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
wxNodalDGgeometry2D<REAL>::Normals2d(int k, REAL norms[])
{
    int f;

    REAL x1 = _xcoord[k][_Fmask[0*_NpF]], y1 = _ycoord[k][_Fmask[0*_NpF]];
    REAL x2 = _xcoord[k][_Fmask[1*_NpF]], y2 = _ycoord[k][_Fmask[1*_NpF]];
    REAL x3 = _xcoord[k][_Fmask[2*_NpF]], y3 = _ycoord[k][_Fmask[2*_NpF]];

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
wxNodalDGgeometry2D<REAL>::GeometricFactors2d(int k, REAL geom[])
{
//    int test = _Fmask[0*_NpF];

    REAL x1 = _xcoord[k][_Fmask[0*_NpF]], y1 =  _ycoord[k][_Fmask[0*_NpF]];
    REAL x2 = _xcoord[k][_Fmask[1*_NpF]], y2 =  _ycoord[k][_Fmask[1*_NpF]];
    REAL x3 = _xcoord[k][_Fmask[2*_NpF]], y3 =  _ycoord[k][_Fmask[2*_NpF]];

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
    geom[0] =  dyds/J;
    geom[1] = -dydr/J;
    geom[2] = -dxds/J;
    geom[3] =  dxdr/J;
    geom[4] =  J;
}

template <typename REAL>
void
wxNodalDGgeometry2D<REAL>::multiplyBYinverseMassMatrix(REAL* input, REAL* output)
{
    for(unsigned kk=0; kk<_NpE*_meqn; kk++)
        output[kk] = 0.0;

    for(unsigned kx=0; kx<_NpE; kx++)
        for(unsigned ky=0; ky<_NpE; ky++)
            for(unsigned kz=0; kz<_meqn; kz++)
                output[kx*_meqn+kz] += _VVT[kx*_NpE+ky]*input[ky*_meqn+kz];

}

// instantiations
template class wxNodalDGgeometry2D<float>;
template class wxNodalDGgeometry2D<double>;
