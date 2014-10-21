#ifndef WXPDGGEOMETRY_H
#define WXPDGGEOMETRY_H

// WarpX lib includes
#include <wxindexer.h>
#include <petscdmplex.h>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

// std includes
#include <vector>

template <typename REAL>
class WxpDGGeometry
{
  public:
/**
 * Create weights, absciae and Legendre polynomials of given spatial
 * order
 *
 * @param meqn Number of equations
 * @param spatialOrder Spatial order
 */
    WxpDGGeometry(DM dm, unsigned meqn, unsigned Spor);

/** Dtor */
    virtual ~WxpDGGeometry();


  private:
/**
 * Find which faces in one element connect to what
 * face in the neighboring element.
 */
    void FacePair2d(DM dm);

/**
 * Calculate the node coordinates for each element.
 */
    void CalculateNodeCoordinates2d(DM dm);

    typedef boost::numeric::ublas::matrix<REAL> matrix_REAL;
    typedef boost::numeric::ublas::matrix<REAL> matrix_INT;

/** No of equations */
    unsigned _meqn;
/** Spatial order */
    unsigned _SpOr;
/** Data management */
    DM _dm;

unsigned _NpE; // number of points/nodes per element
unsigned _NpF; // number of points/nodes per face
unsigned _NfE;  // number of faces per element
REAL *_r, *_s, *_t;  // (r,s,t) coordinates of reference nodes
matrix_REAL _Ds, _Dr, _LIFT; // element matrices
matrix_INT _Fmask;
REAL **p_Dr, **p_Ds, **p_LIFT;
int **p_Fmask;
unsigned **_EtoV; // element to vertecies connectivity matrix
unsigned **_EToE; /* element to neighbor element (elements numbered by their proc) */
unsigned **_EToF; /* element to neighbor face    (element local number 0,1,2) */
unsigned _Klocal; // number of elements in this processor
unsigned _Vlocal; // number of nodes in this processor
REAL **_xcoord; // node x-coordinates
REAL **_ycoord; // node y-coordinates
REAL **_zcoord; // node z-coordinates

};

#endif // WXPDGGEOMETRY_H
