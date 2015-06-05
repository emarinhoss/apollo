#ifndef WXMAXWELLTRANVERSEMAGNETICBC_H
#define WXMAXWELLTRANVERSEMAGNETICBC_H

// WarpX subsolver includes
#include <wxgridbc.h>

/**
 * Base class which applies boundary conditions to fields living on a
 * grid box.
 */
template <typename REAL>
class WxMaxwellTransverseMagnetic : public WxGridBC<REAL>
{
  public:

/**
 * Construct a new grid-bc object
 */
    WxMaxwellTransverseMagnetic()
       : WxGridBC<REAL>("maxwellTransverseMagneticBC") {
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
    void applyBC(REAL *xc, REAL *nx, REAL *q, REAL *qaux, REAL *qBC);

  private:

};

#endif // WXMAXWELLTRANVERSEMAGNETICBC_H
