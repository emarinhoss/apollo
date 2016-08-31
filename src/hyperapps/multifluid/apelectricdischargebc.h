#ifndef APELECTRICDISCHARGEBC_H
#define APELECTRICDISCHARGEBC_H

// WarpX subsolver includes
#include <wxgridbc.h>

/**
 * Base class which applies boundary conditions to fields living on a
 * grid box.
 */
template <typename REAL>
class ApElectricDischargeBC : public WxGridBC<REAL>
{
  public:

/**
 * Construct a new grid-bc object
 */
    ApElectricDischargeBC()
       : WxGridBC<REAL>("electricDischargeBC") {
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
    REAL _Tinfty, _gammaE, _kv, _Efield, _q, _me, _mi, _mn, _gamma, _kB;



};

#endif // APELECTRICDISCHARGEBC_H
