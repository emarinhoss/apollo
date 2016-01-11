// WarpX subsolver includes
#include "wxeulerwallbc.h"

template <typename REAL>
void
WxEulerWallBC<REAL>::setup(const WxCryptSet& wxc, DM dm)
{
  // call base class setup
  WxGridBC<REAL>::setup(wxc, dm);

}

template <typename REAL>
void
WxEulerWallBC<REAL>::applyBC(REAL *xc, REAL *nx, REAL *q, REAL *qaux, REAL *AreaInts, REAL *qBC)
{
    qBC[0] = q[0];
    qBC[1] = q[1]-2.*(nx[0]*q[1]+nx[1]*q[2])*nx[0];
    qBC[2] = q[2]-2.*(nx[0]*q[1]+nx[1]*q[2])*nx[1];
    qBC[3] = 0.0;
    qBC[4] = q[4];
}

// instantiations
template class WxEulerWallBC<float>;
template class WxEulerWallBC<double>;
