#ifndef APMAXWELLRMFANTENNABC_H
#define APMAXWELLRMFANTENNABC_H

// WarpX subsolver includes
#include <wxgridbc.h>

/**
 * Base class which applies boundary conditions to fields living on a
 * grid box.
 */
template <typename REAL>
class WxMaxwellRMFAntennaBC : public WxGridBC<REAL>
{
  public:

/**
 * Construct a new grid-bc object
 */
    WxMaxwellRMFAntennaBC()
       : WxGridBC<REAL>("maxwellRMFAntennaBC") {
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
    REAL _omega, _B0, _phase, _rise, _pi, _epsilon_r, _mu_r;

};
#endif // APMAXWELLRMFANTENNABC_H
