// WarpX subsolver includes
#include "wxeulerzerogradientbc.h"

template <typename REAL>
void
WxEulerZeroGradientBC<REAL>::setup(const WxCryptSet& wxc, DM dm)
{
  // call base class setup
  WxGridBC<REAL>::setup(wxc, dm);

}

template <typename REAL>
void
WxEulerZeroGradientBC<REAL>::applyBC(REAL *xc, REAL *nx, REAL *q, REAL *qaux, REAL *AreaInts, REAL *qBC)
{
    qBC[0] = q[0];
    qBC[1] = q[1];
    qBC[2] = q[2];
    qBC[3] = q[3];
    qBC[4] = q[4];
}

// instantiations
//template class WxEulerZeroGradientBC<float>;
template class WxEulerZeroGradientBC<double>;
