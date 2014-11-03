#ifndef WXPDGGEOMETRY_H
#define WXPDGGEOMETRY_H

// WarpX lib includes
#include <wxindexer.h>
#include <petscdmplex.h>
//#include <boost/numeric/ublas/vector.hpp>
//#include <boost/numeric/ublas/matrix.hpp>

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

/** Destroctor */
    virtual ~WxpDGGeometry();

/** return the number of nodes per element */
    unsigned NpElem(){
        return _NpE;
    }

/** return the x-coordinate for element K node N */
    REAL Xcoordinate(unsigned K, unsigned N){
        return _xcoord[K][N];
    }

/** return the x-coordinate for element K node N */
    REAL Ycoordinate(unsigned K, unsigned N){
        return _ycoord[K][N];
    }

/** calculate the  geometric factors for element k*/
    void GeometricFactors2d(int k, REAL *drdx, REAL *dsdx, REAL *drdy, REAL *dsdy, REAL *J);

/** calculate the face normals for a given element k*/
    void Normals2d(int k, REAL *nx, REAL *ny, REAL *sJ);

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

//    typedef boost::numeric::ublas::matrix<REAL> matrix_REAL;
//    typedef boost::numeric::ublas::matrix<REAL> matrix_INT;

/** No of equations */
    unsigned _meqn;
/** Spatial order */
    unsigned _SpOr;
/** Data management */
    DM _dm;

    int _NpE; // number of points/nodes per element
    int _NpF; // number of points/nodes per face
    int _NfE;  // number of faces per element
    REAL *_r, *_s, *_t;  // (r,s,t) coordinates of reference nodes
    REAL *_Ds, *_Dr,*_LIFT; // element matrices
    int *_Fmask;
    int **_EtoV; // element to verticies connectivity matrix
    int **_EToE; /* element to neighbor element (elements numbered by their proc) */
    int **_EToF; /* element to neighbor face    (element local number 0,1,2) */
    int _Klocal; // number of elements in this processor
    int _Vlocal; // number of nodes in this processor
    REAL **_xcoord; // node x-coordinates
    REAL **_ycoord; // node y-coordinates
    REAL **_zcoord; // node z-coordinates

};

#endif // WXPDGGEOMETRY_H
