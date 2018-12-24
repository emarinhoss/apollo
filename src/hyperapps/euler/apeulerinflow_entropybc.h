#ifndef APEULERINFLOW_ENTROPYBC_H
#define APEULERINFLOW_ENTROPYBC_H


// WarpX subsolver includes
#include <wxgridbc.h>

/**
 * Base class which applies boundary conditions to fields living on a
 * grid box.
 */
template <typename REAL>
class WxEulerInflowEntropyBC : public WxGridBC<REAL>
{
  public:

/**
 * Construct a new grid-bc object
 */
    WxEulerInflowEntropyBC()
       : WxGridBC<REAL>("eulerInflowEntropyBC") {
    }

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
    void applyBC(REAL *xc, REAL *nx, REAL *q, REAL *qaux, REAL *AreaInts, REAL *qBC);

  private:
    REAL _Mach, _rho, _angle, _pres, _gamma;

};

#endif // APEULERINFLOW_ENTROPYBC_H
