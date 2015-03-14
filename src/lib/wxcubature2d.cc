#include "wxcubature2d.h"
// WarpX lib includes
#include <wxmath.h>
#include <wxlogger.h>
#include <wxlogstream.h>
#include "wxcubaturedata2d.h"
#include "wxNodalDGMatrices.h"

template <typename REAL>
WxCubature2d<REAL>::WxCubature2d(DM dm, unsigned meqn, unsigned polOrd, Mat invV)
    : _meqn(meqn), _polyOrd(polOrd), _dm(dm), _iV(invV)
{

    WxLogger *log = WxLogger::get("apollo-root.console");
    WxLogStream infStrm = log->getInfoStream();

    _NPE = (_polyOrd+1)*(_polyOrd+2)/2;
    // Cubature points
    _cubOrd = (int)floor(3.0*(_polyOrd+1));
    // Gaussian quadrature points
    _gQuad = (int)2*ceil((_polyOrd+1)/2.0);

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
            infStrm << "## Error: Invalid order for 2D cubature.  Max polynimoal order should be 8. ##" << std::endl;
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

    // evaluate generalized Vandermonde of Lagrange interpolation functions at cubature nodes
    MatCreateSeqDense(PETSC_COMM_SELF,_pts,_NPE,PETSC_NULL,&Vout);
    Vandermonde2D(_polyOrd,_pts,_r,_s,Vout);
    MatMatMult(Vout,_iV,MAT_INITIAL_MATRIX, PETSC_DEFAULT,&V);
//    MatView(V,PETSC_VIEWER_STDOUT_WORLD);
    // and its transpose
    MatTranspose(V,MAT_INITIAL_MATRIX,&VT);
//    MatView(VT,PETSC_VIEWER_STDOUT_WORLD);
    infStrm << "** done -- Creating Cubature Vandermonde Matrix **" << std::endl;

    // evaluate local derivatives of Lagrange interpolation at cubature points
    MatCreateSeqDense(PETSC_COMM_SELF,_pts,_NPE,PETSC_NULL,&Dr);
    MatCreateSeqDense(PETSC_COMM_SELF,_pts,_NPE,PETSC_NULL,&Ds);
    DifferentiationMatrices2D(_polyOrd,_pts,_r,_s,_iV,&Dr,&Ds);
//    MatView(Dr,PETSC_VIEWER_STDOUT_WORLD);

    // and their transpose
    MatrixTranspose(Dr,&DrT);
    MatrixTranspose(Ds,&DsT);
    infStrm << "** done -- Creating Cubature differentiation Matrices **" << std::endl;

    /** Evaluate data for the Gauss-Legendre quadrature points needed
     * at the element boundaries */
    infStrm << "** Number of Gaussian points per Edge : " << _gQuad << std::endl;
    _gw = alloc_1d<REAL>(_gQuad);
    _gz = alloc_1d<REAL>(_gQuad);
    // weights and abscissa for Gaussian quadrature
    gauleg<REAL>(-1.0, 1.0, _gz-1, _gw-1, _polyOrd);

    int facer[3*_gQuad], faces[3*_gQuad];
    for(unsigned k=0; k<_gQuad; k++)
    {
        facer[0*_gQuad+k] =  _gz[k]; faces[0*_gQuad+k] = -1.0;
        facer[1*_gQuad+k] = -_gz[k]; faces[1*_gQuad+k] =  _gz[k];
        facer[2*_gQuad+k] = -1.0;    faces[2*_gQuad+k] = -_gz[k];
    }

    Mat V;
    MatCreateSeqDense(PETSC_COMM_SELF,3*_gQuad,_NPE,PETSC_NULL,&V);

    Vandermonde2D(_polyOrd,3*_gQuad,facer,faces,V);

    // interpolate at the elemenet edges
    MatMatMult(V,_iV,MAT_INITIAL_MATRIX, PETSC_DEFAULT,&interp);
    // inverse interpolation
    MatrixTranspose(interp,&interpT);
    infStrm << "** done -- Creating Surface interpolation matrices **" << std::endl;
}

template <typename REAL>
WxCubature2d<REAL>::~WxCubature2d()
{
    delete [] _r;
    delete [] _s;
    delete [] _w;
    delete [] _gz;
    delete [] _gw;
    MatDestroy(&V);
    MatDestroy(&Dr);
    MatDestroy(&Ds);
    MatDestroy(&DrT);
    MatDestroy(&DsT);
    MatDestroy(&VT);
    MatDestroy(&Vout);
    MatDestroy(&interp);
    MatDestroy(&interpT);
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
WxCubature2d<REAL>::interpolatedTOCubatures(Vec input, Vec output)
{
    Mat intExp;
    MatCreateSeqAIJ(PETSC_COMM_SELF,_pts*_meqn,_NPE*_meqn,_NPE,PETSC_NULL,&intExp);
    numEqnMatExpand(_pts,V,&intExp);
    MatrixVectorMult(intExp,input,&output);
}

template <typename REAL>
void
WxCubature2d<REAL>::evaluatedVolumeIntegrals(Vec xcoords, Vec ycoords, Vec Fflux, Vec Gflux, Vec *VolInt)
{
    // calculate geometric factors
    Vec rx, sx, ry, sy, J;

    VecCreateSeq(PETSC_COMM_SELF,_pts,&rx);
    VecDuplicate(rx,&sx);
    VecDuplicate(rx,&ry);
    VecDuplicate(rx,&sy);
    VecDuplicate(rx,&J);
    geometricFactors2D(xcoords,ycoords,&rx,&sx,&ry,&sy,&J);

    numEqnVecExpand(&sx);
    numEqnVecExpand(&ry);
    numEqnVecExpand(&sy);
    numEqnVecExpand(&J);

    // cubature weights
    Vec weights;
    PetscScalar *xx;
    VecCreateSeq(PETSC_COMM_SELF,_pts*_meqn,&weights);
    VecGetArray(weights,&xx);
    for(unsigned kk=0; kk<_pts;kk++)
        for(unsigned bb=0; bb<_meqn; bb++)
            xx[kk*_meqn+bb] = _w[kk];
    VecRestoreArray(weights,&xx);

    // Evaluate derivatives
    Vec ddr; VecDuplicate(weights,&ddr);
    evalDerivatives(weights,rx,ry,Fflux,Gflux,DrT,&ddr);
    evalDerivatives(weights,sx,sy,Fflux,Gflux,DsT,VolInt);
    VecAXPY(*VolInt,1.0,ddr);

//    VecDestroy(&sx);
//    VecDestroy(&ry);
//    VecDestroy(&rx);
//    VecDestroy(&sy);

}

template <typename REAL>
void
WxCubature2d<REAL>::evalDerivatives(Vec W, Vec XX, Vec YY, Vec F, Vec G, Mat DD, Vec *DX)
{
    Vec v1, v2;
    VecDuplicate(XX,&v1);
    VecDuplicate(XX,&v2);

    VecPointwiseMult(v1,XX,F);
    VecPointwiseMult(v2,YY,G);
    VecAXPY(v1,1.0,v2);
    VecPointwiseMult(v2,v1,W);
    MatrixVectorMult(DD,v2,DX);

    // destroy
    VecDestroy(&v1);
    VecDestroy(&v2);
}

template <typename REAL>
void
WxCubature2d<REAL>::geometricFactors2D(Vec xcoords, Vec ycoords, Vec *rx, Vec *sx, Vec *ry, Vec *sy, Vec *J)
{
    Vec xr, xs, yr, ys, v1;

    VecCreateSeq(PETSC_COMM_SELF,_pts,&xr);
    VecDuplicate(xr,&xs);
    VecDuplicate(xr,&yr);
    VecDuplicate(xr,&ys);
    VecDuplicate(xr,&v1);

    MatrixVectorMult(Dr,xcoords,&xr);
    MatrixVectorMult(Dr,ycoords,&yr);
    MatrixVectorMult(Ds,xcoords,&xs);
    MatrixVectorMult(Ds,ycoords,&ys);

    VecPointwiseMult(v1,xr,ys);
    VecPointwiseMult(*J,xs,yr);
    VecAYPX(*J,-1.0,v1);

    VecPointwiseDivide(*rx,ys,*J);
    VecPointwiseDivide(*sx,yr,*J); VecScale (*sx, -1.0);
    VecPointwiseDivide(*ry,xs,*J); VecScale (*ry, -1.0);
    VecPointwiseDivide(*sy,xr,*J);

    // destroy
    VecDestroy(&v1);
    VecDestroy(&ys);
    VecDestroy(&yr);
    VecDestroy(&xs);
    VecDestroy(&xr);
}

template <typename REAL>
void
WxCubature2d<REAL>::nodesTOSurfaceGaussians(Vec input, Vec *outPut)
{
    MatrixVectorMult(interp,input,outPut);
}

template <typename REAL>
void
WxCubature2d<REAL>::calculateSurfaceIntegral(Vec numFlux, Vec surfInt)
{
    // cubature weights
    Vec weights, v1;
    PetscScalar *xx;
    VecCreateSeq(PETSC_COMM_SELF,3*_gQuad*_meqn,&weights);
    VecDuplicate(weights,&v1);
    VecGetArray(weights,&xx);
    for(unsigned face=0; face<3; face++)
        for(unsigned kk=0; kk<_gQuad;kk++)
            for(unsigned bb=0; bb<_meqn; bb++)
                xx[(face*_gQuad+kk)*_meqn+bb] = _gw[kk];
    VecRestoreArray(weights,&xx);

    // inverse interpolation matrix
    Mat invInterp;
    MatCreateSeqDense(PETSC_COMM_SELF,_NPE*_meqn,3*_gQuad*_meqn,PETSC_NULL,&invInterp);
    numEqnMatExpand(_NPE,interpT,&invInterp);

    // evaluate integral
    VecPointwiseMult(v1,weights,numFlux);
    MatrixVectorMult(invInterp,v1,&surfInt);

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
WxCubature2d<REAL>::MatrixVectorMult(Mat A, Vec x, Vec *y)
{
    PetscInt mcols, mrows, vrows;
    MatGetSize(A, &mrows, &mcols);

    VecGetSize(x,&vrows);

    WxLogger *log = WxLogger::get("apollo-root.console");
    WxLogStream debugStrm = log->getDebugStream();

    if(mcols!=vrows){
        debugStrm << "** Matrix-Vector Multiplication failed:  " <<
                     mcols << " != " << vrows << std::endl;
        exit(1);
    }

    PetscInt ncols;
    const PetscInt    *cols;
    const PetscScalar *vals;
    PetscScalar *xx, *yy;

    VecGetArray(x,&xx);
    VecGetArray(*y,&yy);
    for(unsigned kk=0; kk<mrows; kk++)
    {
        yy[kk] = 0.0;
        MatGetRow(A,kk,&ncols,&cols,&vals);
        for(unsigned kx=0; kx<ncols; kx++)
            yy[kk] += vals[kx]*xx[kx];
    }
    VecRestoreArray(x,&xx);
    VecRestoreArray(*y,&yy);

}

// instantiations
template class WxCubature2d<float>;
template class WxCubature2d<double>;
