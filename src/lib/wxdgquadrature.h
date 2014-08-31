#ifndef WXDGQUADRATURE_H
#define WXDGQUADRATURE_H

// WarpX lib includes
#include <wxindexer.h>

// std includes
#include <vector>

template <typename REAL>
class WxDGQuadrature
{
  public:
/**
 * Create weights, absciae and Legendre polynomials of given spatial
 * order
 *
 * @param meqn Number of equations
 * @param spatialOrder Spatial order
 */
    WxDGQuadrature(unsigned meqn, unsigned spatialOrder, bool useGaussian=true);

/** Dtor */
    virtual ~WxDGQuadrature();

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
 * Return legendre polynomial at quadrature point
 *
 * @param lm Order of legendre polynomials
 * @param cc Index of quadrature point
 * @return Legenrdre polynomial of order lm at quadrature point cc
 */
    REAL legpol(unsigned lm, unsigned cc) const {
      return _legpol[lm][cc];
    }

/**
 * Return derivative of legendre polynomial at quadrature point
 *
 * @param lm Order of legendre polynomials
 * @param cc Index of quadrature point
 * @return Derivative of Legenrdre polynomial of order lm at quadrature point cc
 */
    REAL dlegpol(unsigned lm, unsigned cc) const {
      return _dlegpol[lm][cc];
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
 * Evaluate expansion along lower edge
 *
 * @param coeffs array of coefficients
 * @param val computed value
 */
    void evalExpansionLower1D(const REAL coeffs[], REAL val[]);

/**
 * Evaluate expansion along upper edge
 */
    void evalExpansionUpper1D(const REAL coeffs[], REAL val[]);

/**
 * Evaluate 1D expansion
 */
    void evalExpansion1D(int lm,
      const REAL coeffs[], REAL res[]);

/**
 * Evaluate expansion along lower edge
 *
 * @param dir direction edge is perpendicular to
 * @param l basis function index
 * @param coeffs array of coefficients
 * @param val computed value
 */
    void evalExpansionLower2D(unsigned dir, unsigned l,
      const REAL coeffs[], REAL val[]);

/**
 * Evaluate expansion along upper edge
 */
    void evalExpansionUpper2D(unsigned dir, unsigned l,
      const REAL coeffs[], REAL val[]);

/**
 * Evaluate 2D expansion
 */
    void evalExpansion2D(int lm, int ln,
      const REAL coeffs[], REAL res[]);

/**
 * Evaluate expansion along lower edge
 *
 * @param dir direction edge is perpendicular to
 * @param l,s basis function index
 * @param coeffs array of coefficients
 * @param val computed value
 */
    void evalExpansionLower3D(unsigned dir, unsigned l,
      unsigned s, const REAL coeffs[], REAL val[]);

/**
 * Evaluate expansion along upper edge
 */
    void evalExpansionUpper3D(unsigned dir, unsigned l,
      unsigned s, const REAL coeffs[], REAL val[]);

/**
 * Evaluate 3D expansion
 */
    void evalExpansion3D(int lm, int ln, int lp,
      const REAL coeffs[], REAL res[]);

  private:
/** No of equations */
    unsigned _meqn;
/** Spatial order, number of coefficients per element, number of quadrature points */
    unsigned _spatialOrder, _numCoeffPerElem, _numQuadPs;
/** Weigths and abscissa for Gaussian quadrature */
    REAL *_w, *_x;
/** Legendre polynomial at quadrature points */
    REAL **_legpol;
/** Derivative of Legendre polynomial at quadrature points */
    REAL **_dlegpol;
/** Normalization constants for basis function */
    REAL *_Cconst, **_Cconst2D, ***_Cconst3D;
/** Indexer for 1D, 2D and 3D expansions */
    WxIndexer<> _ind1D, _ind2D, _ind3D;
};


#endif // WXDGQUADRATURE_H
