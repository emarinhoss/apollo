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

/** Interpolate nodal values into cubature points
 *  @param input  [in]  - conserved variables at nodal points
 *  @param output [out] - conserved variables at cubature point
 */
    void interpolatedTOCubatures(REAL *input, REAL *output);

/** Return the total number of cubature points used. */
    int numCubaturePoints(){
        return _pts;
    }

/** Return the number of Gaussian points per edge/face. */
    int numGaussianPoints(){
        return _gQuad;
    }

/** Evaluate the volume integrals using cubature integration */
    void evaluatedVolumeIntegrals(REAL *xcoords, REAL *ycoords, REAL *Fflux, REAL *Gflux, REAL *VolInt);

/** Evaluate the surface integral */
    void calculateSurfaceIntegral(REAL *numFlux, REAL *surfInt);

/** Interpolate the nodal values to the Gaussian points at all 3 edges of the element.
 *  The solution is stack in order of the surface number */
    void nodesTOSurfaceGaussians(REAL *input, REAL *outPut);

/** Checks if any of the values in the array are NAN's
 * @param n - array size
 * @param y - array
 * @param msg - location message to help with debugging
 */
    void checkNAN(int n, REAL *y , std::string msg);

  private:
/** create amatrix that takes into account the number of equations in the system */
    void numEqnMatExpand(int N, Mat A, Mat *B);

/** create amatrix that takes into account the number of equations in the system */
    void numEqnVecExpand(Vec *A);

/** calculate the geometric factor */
    void geometricFactors2D(REAL *xcoords, REAL *ycoords, REAL *rx, REAL *sx, REAL *ry, REAL *sy, REAL *J);

/** Calculate the matrix transpose of a given matrix */
    void MatrixTranspose(Mat A, Mat *A_trans);

/** Transfer all the Matrices into row major arrays */
    void petscMatTOArray(Mat A, REAL *array);

/** Matrix Vector multiplication. Petsc has the MatMult function, however
 * this functions is not well behaved for matrices that are not square
 *
 * Ax = y
 */
    void MatrixVectorMult(int rows, int cols, int meqn, REAL *A, REAL *x, REAL *y);

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
    REAL *_W, *_iV; //
    REAL *_V, *_Dr, *_Ds, *_VT, *_DrT, *_DsT; // Matrices evaluated at the cubature points
    Mat inverseV;

/** Gaussian data */
    REAL *_gz, *_gw;
    REAL *_interp, *_interpT;

};

#endif // WXCUBATURE2D_H
