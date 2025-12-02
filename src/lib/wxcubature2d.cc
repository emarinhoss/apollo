#include "wxcubature2d.h"
// WarpX lib includes
#include <wxmath.h>
#include <wxlogger.h>
#include <wxlogstream.h>
#include "wxcubaturedata2d.h"
#include "wxNodalDGMatrices.h"

// Include BLAS/LAPACK for optimized matrix operations
#ifdef USE_BLAS
extern "C" {
    #include <cblas.h>
}
#endif

template <typename REAL>
WxCubature2d<REAL>::WxCubature2d(DM dm, unsigned meqn, unsigned polOrd, Mat invV)
    : _meqn(meqn), _polyOrd(polOrd), _dm(dm), inverseV(invV)
{

    WxLogger *log = WxLogger::get("apollo-root.console");
    WxLogStream infStrm = log->getInfoStream();
    WxLogStream debStrm = log->getDebugStream();

    _NPE = (_polyOrd+1)*(_polyOrd+2)/2;
    // Cubature points
    _cubOrd = (int)floor(3.0*(_polyOrd+1));
    // Gaussian quadrature points
    //_gQuad = (int)ceil((2*_polyOrd+1)/2.0);
    _gQuad = floor((_polyOrd+1)*2);

    // get cubature points and weigths
    int cubPoints[28] = { 1, 3, 6, 6, 7,12,15, 16, 19, 25, 28, 36, 40, 46,
                         54,58,66,73,82,85,93,100,106,118,126,138,145,225};
    _pts = cubPoints[_cubOrd-1];
    REAL CT[3*_pts];
    switch (_cubOrd)
    {
        case  1: cub2D_1(_pts,CT); break;
        case  2: cub2D_2(_pts,CT); break;
        case  3: cub2D_3(_pts,CT); break;
        case  4: cub2D_4(_pts,CT); break;
        case  5: cub2D_5(_pts,CT); break;
        case  6: cub2D_6(_pts,CT); break;
        case  7: cub2D_7(_pts,CT); break;
        case  8: cub2D_8(_pts,CT); break;
        case  9: cub2D_9(_pts,CT); break;
        case 10: cub2D_10(_pts,CT); break;
        case 11: cub2D_11(_pts,CT); break;
        case 12: cub2D_12(_pts,CT); break;
        case 13: cub2D_13(_pts,CT); break;
        case 14: cub2D_14(_pts,CT); break;
        case 15: cub2D_15(_pts,CT); break;
        case 16: cub2D_16(_pts,CT); break;
        case 17: cub2D_17(_pts,CT); break;
        case 18: cub2D_18(_pts,CT); break;
        case 19: cub2D_19(_pts,CT); break;
        case 20: cub2D_20(_pts,CT); break;
        case 21: cub2D_21(_pts,CT); break;
        case 22: cub2D_22(_pts,CT); break;
        case 23: cub2D_23(_pts,CT); break;
        case 24: cub2D_24(_pts,CT); break;
        case 25: cub2D_25(_pts,CT); break;
        case 26: cub2D_26(_pts,CT); break;
        case 27: cub2D_27(_pts,CT); break;
        case 28: cub2D_28(_pts,CT); break;
        default:
            infStrm << "## Error: Invalid order for 2D cubature.  Max polynomial order should be 8. ##" << std::endl;
            break;
    }

    infStrm << "** Number of cubature points used : " << _pts << std::endl;

    _r = alloc_1d<REAL>(_pts);
    _s = alloc_1d<REAL>(_pts);
    _w = alloc_1d<REAL>(_pts);
    for(unsigned k=0; k<_pts; k++)
    {
        _r[k] = CT[3*k+0];
        _s[k] = CT[3*k+1];
        _w[k] = CT[3*k+2];
    }

    Mat Vout, V, gV, VT, Dr, Ds, DrT, DsT, interp, interpT, cubMass;

    // evaluate generalized Vandermonde of Lagrange interpolation functions at cubature nodes
    MatCreateSeqDense(PETSC_COMM_SELF,_pts,_NPE,PETSC_NULL,&Vout);
    Vandermonde2D(_polyOrd,_pts,_r,_s,&Vout);
    MatMatMult(Vout,inverseV,MAT_INITIAL_MATRIX, PETSC_DEFAULT,&V);
//    MatView(V,PETSC_VIEWER_STDOUT_WORLD);
    // and its transpose
    MatTranspose(V,MAT_INITIAL_MATRIX,&VT);
    debStrm << "** done -- Creating Cubature Vandermonde Matrix. **" << std::endl;

    // evaluate local derivatives of Lagrange interpolation at cubature points
    MatCreateSeqDense(PETSC_COMM_SELF,_pts,_NPE,PETSC_NULL,&Dr);
    MatCreateSeqDense(PETSC_COMM_SELF,_pts,_NPE,PETSC_NULL,&Ds);
    MatCreateSeqDense(PETSC_COMM_SELF,_NPE,_pts,PETSC_NULL,&DrT);
    MatCreateSeqDense(PETSC_COMM_SELF,_NPE,_pts,PETSC_NULL,&DsT);
    DifferentiationMatrices2D(_polyOrd,_pts,_r,_s,inverseV,&Dr,&Ds);
//    MatView(Dr,PETSC_VIEWER_STDOUT_WORLD);

    // and their transpose
    MatrixTranspose(Dr,&DrT);
//    MatView(Dr,PETSC_VIEWER_STDOUT_WORLD);
    MatrixTranspose(Ds,&DsT);
//    MatView(Ds,PETSC_VIEWER_STDOUT_WORLD);
    debStrm << "** done -- Creating Cubature differentiation Matrices. **" << std::endl;

    // Calculate inverse of mass matrix
    Mat dummy, VVT;
    MatTranspose(Vout,MAT_INITIAL_MATRIX,&dummy);
    MatMatMult(Vout,dummy,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&VVT);
    MatCreateSeqDense(PETSC_COMM_SELF,_pts,_pts,PETSC_NULL,&cubMass);
    this->invertMatrix(VVT,&cubMass);

    /** Evaluate data for the Gauss-Legendre quadrature points needed
     * at the element boundaries */
    infStrm << "** Number of Gaussian points per Edge : " << _gQuad << std::endl;
    _gw = alloc_1d<REAL>(_gQuad);
    _gz = alloc_1d<REAL>(_gQuad);
    // weights and abscissa for Gaussian quadrature
    gauleg<REAL>(-1.0, 1.0, _gz-1, _gw-1, _gQuad);

    REAL facer[3*_gQuad], faces[3*_gQuad];
    for(unsigned k=0; k<_gQuad; k++)
    {
        facer[0*_gQuad+k] =  _gz[k]; faces[0*_gQuad+k] = -1.0;
        facer[1*_gQuad+k] = -_gz[k]; faces[1*_gQuad+k] =  _gz[k];
        facer[2*_gQuad+k] = -1.0;    faces[2*_gQuad+k] = -_gz[k];
    }

    MatCreateSeqDense(PETSC_COMM_SELF,3*_gQuad,_NPE,PETSC_NULL,&gV);
    Vandermonde2D(_polyOrd,3*_gQuad,facer,faces,&gV);

    // interpolate values of the surface Gaussian points from nodal values
    MatMatMult(gV,inverseV,MAT_INITIAL_MATRIX, PETSC_DEFAULT,&interp);
//    MatView(interp,PETSC_VIEWER_STDOUT_WORLD);
    // reverse interpolation, from surface Gaussian points to nodal values
    MatrixTranspose(interp,&interpT);
    debStrm << "** done -- Creating Surface interpolation matrices. **" << std::endl;

    // allocate memory
    _cMass = alloc_1d<REAL>(_pts*_pts);
    _V    = alloc_1d<REAL>(_pts*_NPE);
    _VT   = alloc_1d<REAL>(_pts*_NPE);
    _Dr  = alloc_1d<REAL>(_pts*_NPE);
    _Ds  = alloc_1d<REAL>(_pts*_NPE);
    _DrT = alloc_1d<REAL>(_pts*_NPE);
    _DsT = alloc_1d<REAL>(_pts*_NPE);
    _interp = alloc_1d<REAL>(3*_gQuad*_NPE);
    _interpT = alloc_1d<REAL>(3*_gQuad*_NPE);

    // store the matrices in row major arrays
    petscMatTOArray(cubMass,_cMass);
    petscMatTOArray(V,_V);
    petscMatTOArray(VT,_VT);
    petscMatTOArray(Dr,_Dr);
    petscMatTOArray(Ds,_Ds);
    petscMatTOArray(DrT,_DrT);
    petscMatTOArray(DsT,_DsT);
    petscMatTOArray(interp,_interp);
    petscMatTOArray(interpT,_interpT);

    // Petsc matrices are no longer needed
    MatDestroy(&V);
    MatDestroy(&VT);
    MatDestroy(&gV);
    MatDestroy(&Dr);
    MatDestroy(&Ds);
    MatDestroy(&DrT);
    MatDestroy(&DsT);
    MatDestroy(&interp);
    MatDestroy(&interpT);
    MatDestroy(&cubMass);
}

template <typename REAL>
WxCubature2d<REAL>::~WxCubature2d()
{
    delete [] _r;
    delete [] _s;
    delete [] _w;
    delete [] _gz;
    delete [] _gw;
    delete [] _V;
    delete [] _Dr;
    delete [] _Ds;
    delete [] _DrT;
    delete [] _DsT;
    delete [] _VT;
    delete [] _cMass;
    delete [] _interp;
    delete [] _interpT;
}

template <typename REAL>
void
WxCubature2d<REAL>::invertMatrix(Mat A, Mat *invA)
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
    // MatLUFactor(inpA,is,is,&iluinfo);
    // Calculate inverse
    MatMatSolve(inpA,B,*invA);

    MatDestroy(&inpA);
    MatDestroy(&B);
}

template <typename REAL>
void
WxCubature2d<REAL>::petscMatTOArray(Mat A, REAL *array)
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
        for(int kx=0; kx<ncols; kx++){
            REAL AA = vals[kx];
            array[sk++] = vals[kx];
        }
    }
}

template <typename REAL>
void
WxCubature2d<REAL>::numEqnMatExpand(int N, Mat A, Mat *B)
{
    PetscInt ncols;
    const PetscInt    *cols;
    const PetscScalar *vals;

    for(unsigned kx=0; kx<N; kx++)
    {
        MatGetRow(A,kx,&ncols,&cols,&vals);
        for(unsigned mx=0; mx<ncols; mx++)
            for(unsigned nx=0; nx<_meqn; nx++)
                MatSetValue(*B, kx*_meqn+nx, cols[mx]*_meqn+nx, vals[mx], INSERT_VALUES);
    }

    MatAssemblyBegin(*B, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(*B, MAT_FINAL_ASSEMBLY);
}

template <typename REAL>
void
WxCubature2d<REAL>::numEqnVecExpand(Vec *A)
{
    Vec BB;
    PetscScalar *xx, *yy;
    PetscInt size;
    VecGetSize(*A,&size);

    VecCreateSeq(PETSC_COMM_SELF,size*_meqn,&BB);

    VecGetArray(*A,&xx);
    VecGetArray(BB,&yy);
    for(unsigned aa=0; aa<size; aa++)
        for(unsigned bb=0; bb<_meqn; bb++)
            yy[aa*_meqn+bb] = xx[aa];
    VecRestoreArray(*A,&xx);
    VecRestoreArray(BB,&yy);

    //VecDestroy(A);
    *A = BB;
    VecDestroy(&BB);
}

template <typename REAL>
void
WxCubature2d<REAL>::interpolatedTOCubatures(int num, REAL* input, REAL* output)
{
    MatrixVectorMult(_pts,_NPE,num,_V,input,output);
}

template <typename REAL>
void
WxCubature2d<REAL>::evaluatedVolumeIntegrals(REAL* xcoords, REAL* ycoords, REAL* Fflux, REAL* Gflux, REAL *Src, REAL *VolInt)
{
//    checkNAN(_pts*_meqn, Fflux, "Not-a-number in Fflux.\n");
//    checkNAN(_pts*_meqn, Gflux, "Not-a-number in Gflux.\n");
//    checkNAN(_pts*_meqn, Src, "Not-a-number in Source.\n");


    REAL rx[_pts], sx[_pts], ry[_pts], sy[_pts], J[_pts];
    geometricFactors2D(xcoords,ycoords,rx,sx,ry,sy,J);

    REAL ddx[_NPE*_meqn], ddy[_NPE*_meqn], src[_NPE*_meqn], RR[_pts*_meqn], SS[_pts*_meqn], CC[_pts*_meqn];

    for(unsigned kk=0; kk<_pts; kk++)
        for(unsigned mm=0; mm<_meqn; mm++){
            RR[kk*_meqn+mm] = _w[kk]*J[kk]*(rx[kk]*Fflux[kk*_meqn+mm]+ry[kk]*Gflux[kk*_meqn+mm]);
            SS[kk*_meqn+mm] = _w[kk]*J[kk]*(sx[kk]*Fflux[kk*_meqn+mm]+sy[kk]*Gflux[kk*_meqn+mm]);
            CC[kk*_meqn+mm] = _w[kk]*J[kk]*Src[kk*_meqn+mm];
        }

    MatrixVectorMult(_NPE,_pts,_meqn,_DrT,RR,ddx);
    MatrixVectorMult(_NPE,_pts,_meqn,_DsT,SS,ddy);
    MatrixVectorMult(_NPE,_pts,_meqn,_VT ,CC,src);

    // Add x-dir and y-dir contributions
    for(unsigned kk=0; kk<_NPE*_meqn; kk++)
        VolInt[kk] = ddx[kk] + ddy[kk] + src[kk];

//    checkNAN(_NPE*_meqn, VolInt, "Not-a-number in VolInt.\n");
}

template <typename REAL>
void
WxCubature2d<REAL>::geometricFactors2D(REAL* xcoords, REAL* ycoords, REAL* rx, REAL* sx, REAL* ry, REAL* sy, REAL* J)
{
    REAL xr[_pts], xs[_pts], yr[_pts], ys[_pts];
    MatrixVectorMult(_pts,_NPE,1,_Dr,xcoords,xr);
    MatrixVectorMult(_pts,_NPE,1,_Dr,ycoords,yr);
    MatrixVectorMult(_pts,_NPE,1,_Ds,xcoords,xs);
    MatrixVectorMult(_pts,_NPE,1,_Ds,ycoords,ys);

    // calculate determinant of the Jacobian
    for(unsigned kk=0; kk<_pts; kk++){
        J[kk] = xr[kk]*ys[kk] - xs[kk]*yr[kk];
        if(J[kk]<=0){
            WxLogger *l = WxLogger::get("apollo-root.console");
            WxLogStream errStrm = l->getErrorStream();
            errStrm << "Error: Jacobian determinant is " << J[kk] << ".\n";
//            PetscFinalize();
            exit(1); // abort execution
        }
    }

    for(unsigned kk=0; kk<_pts; kk++)
    {
        rx[kk] =  ys[kk]/J[kk];
        sx[kk] = -yr[kk]/J[kk];
        ry[kk] = -xs[kk]/J[kk];
        sy[kk] =  xr[kk]/J[kk];
    }
}

template <typename REAL>
void
WxCubature2d<REAL>::nodesTOSurfaceGaussians(int len, REAL* input, REAL *outPut)
{
    MatrixVectorMult(3*_gQuad,_NPE,len,_interp,input,outPut);
}

template <typename REAL>
void
WxCubature2d<REAL>::calculateSurfaceIntegral(REAL* numFlux, REAL* surfInt)
{
//    checkNAN(3*_gQuad*_meqn, numFlux, "Not-a-number in numFlux.\n");
    REAL W[3*_gQuad];

    for(unsigned kk=0; kk<3; kk++)
        for(unsigned mm=0; mm<_gQuad; mm++)
            W[kk*_gQuad+mm] = _gw[mm];

    for(unsigned kk=0; kk<3*_gQuad; kk++)
        for(unsigned mm=0; mm<_meqn; mm++)
            numFlux[kk*_meqn+mm] *= W[kk];

    MatrixVectorMult(_NPE,3*_gQuad,_meqn,_interpT,numFlux,surfInt);
//    checkNAN(_NPE*_meqn, surfInt, "Not-a-number in surfInt.\n");
}

template <typename REAL>
void
WxCubature2d<REAL>::discontinuityDetectorIntegral(REAL* surfaceVals, REAL* IntPerFace)
{
    REAL W[3*_gQuad];

    for(unsigned kk=0; kk<_meqn; kk++)
        IntPerFace[kk] = 0.0;

    for(unsigned kk=0; kk<3; kk++)
        for(unsigned mm=0; mm<_gQuad; mm++)
            W[kk*_gQuad+mm] = _gw[mm];

    for(unsigned kk=0; kk<3*_gQuad; kk++)
        for(unsigned mm=0; mm<_meqn; mm++)
            surfaceVals[kk*_meqn+mm] *= W[kk];

    for(unsigned kk=0; kk<3; kk++)
        for(unsigned mm=0; mm<_meqn; mm++)
            for(unsigned quads=0; quads<_gQuad; quads++)
                IntPerFace[mm] += surfaceVals[kk*_gQuad*_meqn+quads*_meqn+mm];
}


template <typename REAL>
void
WxCubature2d<REAL>::MatrixTranspose(Mat A, Mat *A_trans)
{
    PetscInt ccols, rrows;
    MatGetSize(A, &rrows, &ccols);
    MatCreateSeqDense(PETSC_COMM_SELF,ccols,rrows,PETSC_NULL,A_trans);

    PetscInt ncols;
    const PetscInt    *cols;
    const PetscScalar *vals;

    for(unsigned kk=0; kk<rrows; kk++)
    {
        MatGetRow(A,kk,&ncols,&cols,&vals);
        for(unsigned kx=0; kx<ncols; kx++)
            MatSetValue(*A_trans, kx, kk, vals[kx], INSERT_VALUES);
    }

    MatAssemblyBegin(*A_trans, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(*A_trans, MAT_FINAL_ASSEMBLY);
}

template <typename REAL>
void
WxCubature2d<REAL>::MatrixVectorMult(int rows, int cols, int meqn, REAL *A, REAL *x, REAL *y)
{
    // Use optimized BLAS for matrix-vector multiplication when available
    // This computes: Y = A * X where X and Y have multiple columns (meqn)
    // Equivalent to: for each column i: y[:,i] = A * x[:,i]

#ifdef USE_BLAS
    // BLAS dgemm: C = alpha*A*B + beta*C
    // We use it as: Y(rows x meqn) = A(rows x cols) * X(cols x meqn)
    if(sizeof(REAL) == sizeof(double)) {
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                    rows, meqn, cols,
                    1.0, (double*)A, cols, (double*)x, meqn,
                    0.0, (double*)y, meqn);
    } else if(sizeof(REAL) == sizeof(float)) {
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                    rows, meqn, cols,
                    1.0f, (float*)A, cols, (float*)x, meqn,
                    0.0f, (float*)y, meqn);
    }
#else
    // Fallback to original implementation if BLAS not available
    for(int kk=0; kk<rows*meqn; kk++)
        y[kk] = 0.0;

    for(unsigned kx=0; kx<rows; kx++)
        for(unsigned ky=0; ky<cols; ky++)
            for(unsigned kz=0; kz<meqn; kz++)
                y[kx*meqn+kz] += A[kx*cols+ky]*x[ky*meqn+kz];
#endif
}

template <typename REAL>
void
WxCubature2d<REAL>::checkNAN(int n, REAL *y, std::string msg)
{
    for(int kk=0; kk<n; kk++)
    {
        if(y[kk]!=y[kk])
        {
            WxLogger *l = WxLogger::get("apollo-root.console");
            WxLogStream errStrm = l->getInfoStream();
            errStrm << msg ;
//            PetscFinalize();
            exit(1); // abort execution
        }
    }
}

template <typename REAL>
void
WxCubature2d<REAL>::CalculateAreaIntegrals(REAL* xcoords, REAL* ycoords, REAL *Src, REAL *AreaInt)
{
    REAL rx[_pts], sx[_pts], ry[_pts], sy[_pts], J[_pts];
    geometricFactors2D(xcoords,ycoords,rx,sx,ry,sy,J);

    for(unsigned mm=0; mm<_meqn; mm++)
        AreaInt[mm] = 0.0;

    for(unsigned mm=0; mm<_meqn; mm++)
        for(unsigned kk=0; kk<_pts; kk++)
            AreaInt[mm] += _w[kk]*Src[kk*_meqn+mm]*J[kk];
}

// instantiations
//template class WxCubature2d<float>;
template class WxCubature2d<double>;
