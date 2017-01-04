#ifndef APTUALIABADILIMITER_H
#define APTUALIABADILIMITER_H

// WarpX subsolver includes
#include <wxnodaldglimiter.h>
#include <wxhyperboliceqnset.h>
#include <wxgridbc.h>

/**
 * Base class which applies boundary conditions to fields living on a
 * grid box.
 */
template <typename REAL>
class WxTuAliabadiLimiter : public WxNodalDGLimiter<REAL>
{
  public:

/**
 * Construct a new grid-bc object
 */
    WxTuAliabadiLimiter()
       : WxNodalDGLimiter<REAL>("tuAliabadiLimiter") {
    }

/** Destructor */
    virtual ~WxTuAliabadiLimiter();

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
/** Set of hyperbolic equations to solve */
    WxHyperbolicEqnSet<REAL> _eqnSet;
    unsigned _meqn;

    DM _dm;         // data structure/management object
    int _Klocal;    // Total number of element

/** Gradient data is stored */
    // Edge gradients
    REAL *dVdxE1, *dVdyE1, *dVdxE2, *dVdyE2, *dVdxE3, *dVdyE3;

/** Array to modify the directions when doing cartesian v/s radial */
    unsigned _dirs[3];

    // Cell average gradients
    REAL *dVdxC0, *dVdyC0;

    std::vector<std::string> _bcSubSolvers; // boundary conditions

/**
 *  Apply boundary conditions
 */
    void applyBc(int bcNum, REAL *xc, REAL *nx, REAL *q, REAL *qaux, REAL *AreaInts, REAL *qBC);

/**
 * Infinity/Not-a-Number check
 */
    PetscErrorCode isInfinityOrNAN(Vec f, std::string location);


};


#endif // APTUALIABADILIMITER_H
