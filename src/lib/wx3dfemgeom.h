#ifndef WX3DFEMGEOM_H
#define WX3DFEMGEOM_H

// WarpX lib includes
#include <wxindexer.h>
#include <wxgridbox.h>

// std includes
#include <vector>

template <typename REAL>
class Wx3dFEMGeom
{
  public:
/**
 * Create weights, absciae and Legendre polynomials of given spatial
 * order
 *
 * @param meqn Number of equations
 * @param spatialOrder Spatial order
 */
    Wx3dFEMGeom(unsigned meqn, unsigned femSpor, unsigned dgSpor, WxGridBox<REAL> gb);

/** Dtor */
    virtual ~Wx3dFEMGeom();

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
* Evaluate conserved values at a given quadrature point
* @param li - x quadrature index
* @param lj - y quadrature index
* @param lk - z quadrature index
* @param q1 - conserved values at element nodes
* @param res- conserved values at specified quadrature locations
*/
    void evalQvalue3D(int li, int lj, int lk, REAL **q1, REAL *res);


/**
 * Create the x-coordinates of the nodes
 *
 * @param gb the grid
 * @param rhs coordinates matrix
 */
    void getNodeCoordinates(WxGridBox<REAL> gb, REAL **rhs);


/**
 * Create the element connectivity matrix
 *
 * @param spOrd spatial order
 * @param rhs connectivity matrix
 */
    void buildConnectivity(int **rhs);

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
 * Maps the ijk coordinates of the element to the
 * element number
 *
 * @return ijk to element number
 * @return element number to ijk
 */
   void ijktoelementnumber(REAL ***rhs, REAL **rhs2);

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

/**
 * Returns the ijk coordinates of a given element number
 *
 */
   void ijkvalues(unsigned num, unsigned ijk[3]) const{
       ijk[0] = _entoijk[num][0];
       ijk[1] = _entoijk[num][1];
       ijk[2] = _entoijk[num][2];
   }

/**
 * Returns the element number from ijk coordinates
 *
 */
   int elemnumber(int i, int j, int k) const{
       return _ijktoen[i][j][k];
   }

  private:
/** No of equations */
    unsigned _meqn;
/** No of spatial order for both methods */
    unsigned _femSpOr, _dgSpOr;
/** Number of quadrature points and nodes/modes per element */
    unsigned _quadPoints, _nodesPerElem, _modesPerElem;
    unsigned _nxNodesElem, _nyNodesElem, _nzNodesElem;
/** Total number of nodes */
    unsigned _TotNumNodes;
/** Total number of elements */
    unsigned _TotNumElem, _TotNumElemNX, _TotNumElemNY, _TotNumElemNZ;
/** Weigths and abscissa for Gaussian quadrature */
    REAL *_w, *_x;
/** Shape function at quadrature points */
    REAL **_Phi, **_legpol;
/** Derivative wrt to natural coords of Shape function at quadrature points */
    REAL **_dPhi , **_dlegpol;
/** Derivative wrt to global coords of Shape function at quadrature points */
    REAL **_dPhi_x;
/** Indexer for 1D, 2D and 3D expansions */
    WxIndexer<> _ind3D;
/** Coordinate of the nodes and connectivity matrix*/
    REAL **_coords, **_index;
    int **_connect, _ndims;
/** element length */
    REAL _dx;
/** Normalization constants for basis function */
    REAL ***_Cconst3D;
/** IJK to element number mapping and vice-versa */
    int ***_ijktoen, **_entoijk;
};

#endif // WX3DFEMGEOM_H
