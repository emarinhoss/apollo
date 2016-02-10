#ifndef WXHESTHAVENWARBURTONEULERLIMITER_H
#define WXHESTHAVENWARBURTONEULERLIMITER_H

// WarpX subsolver includes
#include <wxnodaldglimiter.h>
#include <wxgridbc.h>

/**
 * Base class which applies boundary conditions to fields living on a
 * grid box.
 */
template <typename REAL>
class WxHestavenWarburtonEulerLimiter : public WxNodalDGLimiter<REAL>
{
  public:

/**
 * Construct a new grid-bc object
 */
    WxHestavenWarburtonEulerLimiter()
       : WxNodalDGLimiter<REAL>("eulerLimiterHW") {
    }

/** Destructor */
    virtual ~WxHestavenWarburtonEulerLimiter();

  protected:

/**
 * Setup subsolver object using supplied cryptset
 *
 * @param wxc Cryptset to use for setting
 */
    void setup(const WxCryptSet& wxc, DM dm);

/**
 * Apply BC to the upper edge along direction 'dir'
 *
 * @param dir direction in which to apply BC
 * @param arr array to which apply BC
 */
    void applyLimiter(wxNodalDGgeometry2D<REAL> *geom, WxCubature2d<REAL> *cub, Vec qk, Vec q_limited);

  private:
    REAL _gamma;    // ratio of specific heats for Euler equantions
    DM _dm;         // data structure/management object
    int _Klocal;    // Total number of element
    REAL _minPres, _minDens; // minimum pressure and density values

/** Gradient data is stored */
    // Edge gradients
    REAL *dVdxE1, *dVdyE1, *dVdxE2, *dVdyE2, *dVdxE3, *dVdyE3;

    // Cell average gradients
    REAL *dVdxC0, *dVdyC0;

    std::vector<std::string> _bcSubSolvers; // boundary conditions

    void computeConservedAndPrimitiveAVEVariables(int kNodes, PetscScalar *qIn, REAL *AVE, REAL *qCons, REAL *qPrim);

/**
 *  Apply boundary conditions
 */
    void applyBc(int bcNum, REAL *xc, REAL *nx, REAL *q, REAL *qaux, REAL *AreaInts, REAL *qBC);

/**
 * Compute the primitive variables from the conserved variables
 */
    void primitiveVariables(REAL *qCons, REAL *qPrim);

};

#endif // WXHESTHAVENWARBURTONEULERLIMITER_H
