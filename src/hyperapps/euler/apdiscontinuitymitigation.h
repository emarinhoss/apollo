#ifndef APDISCONTINUITYMITIGATION_H
#define APDISCONTINUITYMITIGATION_H

// WarpX subsolver includes
#include <wxhyperboliceqnset.h>
#include <wxnodaldglimiter.h>
#include <wxgridbc.h>

/**
 * Base class which applies boundary conditions to fields living on a
 * grid box.
 */
template <typename REAL>
class ApDiscontinuityMitigation : public WxNodalDGLimiter<REAL>
{
  public:

/**
 * Construct a new grid-bc object
 */
    ApDiscontinuityMitigation()
       : WxNodalDGLimiter<REAL>("discontinuityMitigation") {
    }

/** Destructor */
    virtual ~ApDiscontinuityMitigation();

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

DM _dm;         // data structure/management object
int _Klocal;    // Total number of element

std::vector<std::string> _bcSubSolvers; // boundary conditions
/** Array to modify the directions when doing cartesian v/s radial */
  unsigned _dirs[3];

/** Set of hyperbolic equations to solve */
  WxHyperbolicEqnSet<REAL> _eqnSet;

/** Equations */
  unsigned _meqn;

/** Polynomial Order */
  unsigned _polyOrder;

  // jumps, cons. var in left and right of edge i
  REAL *_qM, *_qP, *_surfacesIntegral;

/**
 *  Apply boundary conditions
 */
    void applyBc(int bcNum, REAL *xc, REAL *nx, REAL *q, REAL *qaux, REAL *AreaInts, REAL *qBC);

/**
 * Detect the elements that need to apply limiter/artificial dissipation
 */
    void detectElements(wxNodalDGgeometry2D<REAL> *geom, WxCubature2d<REAL> *cub, Vec qk, int *limElems);

};


#endif // APDISCONTINUITYMITIGATION_H
