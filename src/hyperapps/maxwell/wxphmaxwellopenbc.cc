#include "wxphmaxwellopenbc.h"

template <typename REAL>
void
WxPHMaxwellOpenBC<REAL>::setup(const WxCryptSet& wxc, DM dm)
{
  // call base class setup
  WxGridBC<REAL>::setup(wxc, dm);

}

template <typename REAL>
void
WxPHMaxwellOpenBC<REAL>::applyBC(REAL *xc, REAL *nx, REAL *q, REAL *qaux, REAL *AreaInts, REAL *qBC)
{

    for(unsigned k=0; k<8; k++)
        qBC[k] = q[k];

}

// instantiations
//template class WxPHMaxwellOpenBC<float>;
template class WxPHMaxwellOpenBC<double>;
