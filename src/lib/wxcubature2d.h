#ifndef WXCUBATURE2D_H
#define WXCUBATURE2D_H

// WarpX lib includes
#include <wxindexer.h>
#include <petscdmplex.h>

// std includes
#include <vector>

typedef struct{

} Cub2D;

template <typename REAL>
class WxCubature2d
{
  public:
/**
 * Create weights, absciae and Legendre polynomials of given spatial
 * order
 *
 * @param meqn Number of equations
 * @param spatialOrder Spatial order
 */
    WxCubature2d(DM dm, unsigned meqn, unsigned Spor, Mat invV);

/** Destroctor */
    virtual ~WxCubature2d();

/** Interpolate nodal values into cubature points */
    void interpolatedTOCubatures(Vec input, Vec output);

/** Return the number of cubature points */
    int numCubaturePoints(){
        return _pts;
    }
/** Return the number of Gaussian points per face */
    int numGaussianPoints(){
        return _gQuad;
    }

/** Evaluate the volume integrals using cubature integration */
    void evaluatedVolumeIntegrals(Vec xcoords, Vec ycoords, Vec Fflux, Vec Gflux, Vec *VolInt);

/** Evaluate the surface integral */
    void calculateSurfaceIntegral(Vec numFlux, Vec surfInt);

/** Interpolate the nodal values to the Gaussian points at all 3 edges of the element.
 *  The solution is stack in order of the surface number */
    void nodesTOSurfaceGaussians(Vec input, Vec *outPut);

  private:
/** create amatrix that takes into account the number of equations in the system */
    void numEqnMatExpand(int N, Mat A, Mat *B);

/** create amatrix that takes into account the number of equations in the system */
    void numEqnVecExpand(Vec *A);

/** calculate the geometric factor */
    void geometricFactors2D(Vec xcoords, Vec ycoords, Vec *rx, Vec *sx, Vec *ry, Vec *sy, Vec *J);

/** Evaluate X or Y derivatives at cubature points */
    void evalDerivatives(Vec W, Vec XX, Vec YY, Vec F, Vec G, Mat DD, Vec *DX);

/** Calculate the matrix transpose of a given matrix */
    void MatrixTranspose(Mat A, Mat *A_trans);

/** Matrix Vector multiplication. Petsc has the MatMult function, however
 * this functions is not well behaved for matrices that are not square
 *
 * Ax = b
 */
    void MatrixVectorMult(Mat A, Vec x, Vec *y);

/** No of equations */
    unsigned _meqn;
/** Polynimial Order and number of nodes per element */
    unsigned _polyOrd, _NPE;
/** Cubature Order and number of cubature points*/
    unsigned _cubOrd, _pts;
/** Face Gaussian points */
    unsigned _gQuad;
/** Data management */
    DM _dm;

/** Cubature data */
    REAL *_r, *_s, *_w; // cubature coordinates and weights
    Mat W, _iV, Vout; //
    Mat V, Dr, Ds, VT, DrT, DsT; // Matrices evaluated at the cubature points

/** Gaussian data */
    REAL *_gz, *_gw;
    Mat interp, interpT;

};

#endif // WXCUBATURE2D_H
