#ifndef WXFEMGEOMETRY_H
#define WXFEMGEOMETRY_H

// WarpX lib includes
#include <wxindexer.h>
#include <wxgridbox.h>

// std includes
#include <vector>

template <typename REAL>
class WxFEMGeometry
{
  public:
/**
 * Create weights, absciae and Legendre polynomials of given spatial
 * order
 *
 * @param meqn Number of equations
 * @param spatialOrder Spatial order
 */
    WxFEMGeometry(unsigned ieqn, unsigned meqn, unsigned femSpor, unsigned dgSpor, WxGridBox<REAL> gb, unsigned nodes);

/** Dtor */
    virtual ~WxFEMGeometry();

/**
 * Return quadrature point
 *
 * @param cc Index of quadrature point
 * @return Coordinate of quadrature point in range [-1,1]
 */
    REAL x(unsigned cc) const {
      return _x[cc];
    }

/**
 * Return weight at quadrature point.
 *
 * @param cc Index of quadrature point
 * @return Weight at quadrature
 */
    REAL w(unsigned cc) const {
      return _w[cc];
    }

/**
 * Return legendre/lagrange polynomial at quadrature point
 *
 * @param shape function at node cc
 * @param lm index of quadrature point
 * @return Shape function
 */
    REAL legpol(unsigned cc, unsigned lm) const {
        return _legpol[cc][lm];
    }

    REAL lagpol(unsigned cc, unsigned lm) const {
        return _Phi[cc][lm];
    }

/**
 * Return derivative of legendre/lagrange polynomial at quadrature point
 *
 * @param shape function at node cc
 * @param lm Index of quadrature point
 * @return Derivative of shape function
 */
    REAL dlegpol(unsigned cc, unsigned lm) const {
        return _dlegpol[cc][lm];
    }

    REAL dlagpol(unsigned cc, unsigned lm) const {
        return _dPhi[cc][lm];
    }

/**
 * Return normalization of legendre polynomial in 1D
 *
 * @param lm Order of legendre polynomials
 * @return Normalization of legendre polynomial
 */
    REAL normConst1D(unsigned lm) const {
        return _Cconst[lm];
    }

/**
 * Return normalization of legendre polynomial in 2D
 *
 * @param lm Order of legendre polynomials
 * @return Normalization of legendre polynomial
 */
    REAL normConst2D(unsigned lm, unsigned ln) const {
      return _Cconst2D[lm][ln];
    }

/**
 * Return normalization of legendre polynomial in 3D
 *
 * @param lm Order of legendre polynomials
 * @return Normalization of legendre polynomial
 */
    REAL normConst3D(unsigned lm, unsigned ln, unsigned lp) const {
      return _Cconst3D[lm][ln][lp];
    }

/**
 * Return the element lenth
 *
 * @param shape function at node cc
 * @param lm Index of quadrature point
 * @return Derivative of shape function
 */
    REAL dx() const {
        return _dx;
    }

/**
 * Maps physical coordinate to absicca in [-1,1]
 *
 * @param x physical coordinate
 * @param xcell Cell center coordinate
 * @param dx Cell size
 * @return absicca in [-1,1]
 */
    REAL ETA(REAL x, REAL xcell, REAL dx) const {
      return 2.0*(x-xcell)/dx;
    }

/**
 * Maps absicca in [-1,1] to physical coordinate
 *
 * @param eta Absicca
 * @param xcell Cell center coordinate
 * @param dx Cell size
 * @return physical coordinate
 */
    REAL ATE(REAL eta, REAL xcell, REAL dx) const {
      return 0.5*eta*dx+xcell;
    }

/**
 * Evaluate 1D expansion
 */
    void evalExpansion1Dcg(int lm, REAL **q1, REAL *res);
    void evalExpansion1Ddg(int lm, const REAL coeffs[], REAL res[]);
    void evalExpansionUpper1Ddg(const REAL coeffs[], REAL val[]);
    void evalExpansionLower1Ddg(const REAL coeffs[], REAL val[]);

/**
 * Create the x-coordinates of the nodes
 *
 * @param gb the grid
 * @param rhs coordinates matrix
 */
    void getCoordinate1d(WxGridBox<REAL> gb, REAL **rhs);


/**
 * Create the element connectivity matrix
 *
 * @param spOrd spatial order
 * @param rhs connectivity matrix
 */
    void buildConnectivity1d(int spOrd, int **rhs);

/**
 * Return x-coordinate at quadrature point.
 *
 * @param cc Index of the point
 * @param lm Index of the point
 * @return x-coordinate values at point
 */
    REAL coordinate(unsigned cc, unsigned lm) const {
        return _coords[cc][lm];
    }

/**
 * Return connectivity of element cc and node lm
 *
 * @param element number cc
 * @param node number lm
 * @return connectivity value
 */
   REAL connectivity(unsigned cc, unsigned lm) const {
         return _connect[cc][lm];
    }

/**
 * Return the number of quadrature points
 *
 */
   REAL numQuadPoint() const {
         return _quadPoints;
    }

/**
 * Return connectivity of element cc and node lm
 *
 * @return connectivity value
 */
   REAL nodesPerElem() const {
         return _nodesPerElem;
    }

/**
 * Return connectivity of element cc and node lm
 *
 * @return connectivity value
 */
   REAL modesPerElem() const {
         return _modesPerElem;
    }

  private:
/** No of equations */
    unsigned _meqn, _mauxeqn;
/** No of spatial order for both methods */
    unsigned _femSpOr, _dgSpOr;
/** Number of quadrature points and nodes/modes per element */
    unsigned _quadPoints, _nodesPerElem, _modesPerElem;
/** Total number of nodes */
    unsigned _TotNumNodes;
/** Total number of elements */
    unsigned _TotNumElem;
/** Weigths and abscissa for Gaussian quadrature */
    REAL *_w, *_x;
/** Shape function at quadrature points */
    REAL **_Phi, **_legpol;
/** Derivative wrt to natural coords of Shape function at quadrature points */
    REAL **_dPhi , **_dlegpol;
/** Derivative wrt to global coords of Shape function at quadrature points */
    REAL **_dPhi_x;
/** Indexer for 1D, 2D and 3D expansions */
    WxIndexer<> _ind1D, _ind2D, _ind3D;
/** Coordinate of the nodes and connectivity matrix*/
    REAL **_coords, **_index;
    int **_connect, _ndims;
/** element length */
    REAL _dx;
/** Normalization constants for basis function */
    REAL *_Cconst, **_Cconst2D, ***_Cconst3D;

};

#endif // WXFEMGEOMETRY_H
