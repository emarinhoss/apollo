#ifndef __wxfemquadrature__
#define __wxfemquadrature__

// WarpX lib includes
#include <wxindexer.h>
#include <wxgridbox.h>

// std includes
#include <vector>

template <typename REAL>
class WxFEMQuadrature
{
  public:
/**
 * Create weights, absciae and Legendre polynomials of given spatial
 * order
 *
 * @param meqn Number of equations
 * @param spatialOrder Spatial order
 */
    WxFEMQuadrature(unsigned meqn, unsigned nodesPerElement, WxGridBox<REAL> gb, unsigned qvals);

/** Dtor */
    virtual ~WxFEMQuadrature();

/**
 * Return quadrature point
 *
 * @param cc Index of quadrature point
 * @return Coordinate of quadrature point in range [-1,1]
 */
    REAL x(unsigned cc) const {
      return _x[cc];
    }

    REAL xdg(unsigned cc) const {
      return _dgx[cc];
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
 * Return legendre polynomial at quadrature point
 *
 * @param shape function at node cc
 * @param lm index of quadrature point
 * @return Shape function
 */
    REAL legpol(unsigned cc, unsigned lm) const {
        return _legpol[cc][lm];
    }

/**
 * Return derivative of legendre polynomial at quadrature point
 *
 * @param shape function at node cc
 * @param lm Index of quadrature point
 * @return Derivative of shape function
 */
    REAL dlegpol(unsigned cc, unsigned lm) const {
        return _dlegpol[cc][lm];
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
    void evalExpansion1D(int lm, REAL **q1, REAL *res);
    void evalExpansion1Ddg(int auxMeqn, int lm, REAL **q1, REAL *res);

/**
 * Create the x-coordinates of the nodes
 *
 * @param x-coordinates
 */
    void getCoordinate(WxGridBox<REAL> gb, REAL *rhs);

/**
 * Create the connectivity of the elements
 *
 * @param connectivity matrix
 */
    void buildConnectivity(unsigned spOrd, unsigned n, int **rhs);

/**
 * Return x-coordinate at quadrature point.
 *
 * @param cc Index of the point
 * @return x-coordinate values at point
 */
    REAL coordinate(unsigned cc) const {
      return _coords[cc];
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
 * Number of quadrature points
 *
 * @return number of quadrature points
 */
   REAL numQuadPoint() const {
         return _numQuadPs;
    }

  private:
/** No of equations */
    unsigned _meqn;
/** Number of quadrature points and nodes per element */
    unsigned _numQuadPs, _nodesPerElem;
/** Total number of nodes */
    unsigned _numNodes;
/** Total number of elements */
    unsigned _numElem;
/** Weigths and abscissa for Gaussian quadrature */
    REAL *_w, *_x, *_dgw, *_dgx;
/** Shape function at quadrature points */
    REAL **_legpol, **_dglegpol;
/** Derivative of Shape function at quadrature points */
    REAL **_dlegpol, **_dgdlegpol;
/** Indexer for 1D, 2D and 3D expansions */
    WxIndexer<> _ind1D;
/** Coordinate of the vertices and connectivity matrix*/
    REAL *_coords;
    int **_connect;
/** element length */
    REAL _dx;
/** number of ghost cells */
    int _npad;
};

#endif // __wxfemquadrature__
