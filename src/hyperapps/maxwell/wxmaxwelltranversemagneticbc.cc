// WarpX subsolver includes
#include "wxmaxwelltranversemagneticbc.h"

template <typename REAL>
void
WxMaxwellTransverseMagnetic<REAL>::setup(const WxCryptSet& wxc, DM dm)
{
  // call base class setup
  WxGridBC<REAL>::setup(wxc, dm);

}

template <typename REAL>
void
WxMaxwellTransverseMagnetic<REAL>::applyBC(REAL *xc, REAL *nx, REAL *q, REAL *qaux, REAL *AreaInts, REAL *qBC)
{
    qBC[0] =  q[0];
    qBC[1] =  q[1];
    qBC[2] = -q[2];
    qBC[3] =  q[3];
    qBC[4] =  q[4];
    qBC[5] =  q[5];
    qBC[6] =  q[6];
    qBC[7] =  q[7];
}

// instantiations
//template class WxMaxwellTransverseMagnetic<float>;
template class WxMaxwellTransverseMagnetic<double>;
