#include "apmaxwellrmfantennabc.h"
#include <wxmath.h>

template <typename REAL>
void
WxMaxwellRMFAntennaBC<REAL>::setup(const WxCryptSet& wxc, DM dm)
{
  // call base class setup
  WxGridBC<REAL>::setup(wxc, dm);

  _freq = wxc.template get<REAL>("frequency");
  _B0 = wxc.template get<REAL>("B_rmf");
  _phase = wxc.template get<REAL>("phase");
  _rise = wxc.template get<REAL>("rise_time");

  _pi = 3.141592653589793;

}

template <typename REAL>
void
WxMaxwellRMFAntennaBC<REAL>::applyBC(REAL *xc, REAL *nx, REAL *q, REAL *qaux, REAL *AreaInts, REAL *qBC)
{

    REAL ex = q[0];
    REAL ey = q[1];
    REAL ez = q[2];

    REAL enorm = ex*nx[0] + ey*nx[1];
    REAL etang = ex*nx[1] - ey*nx[0];

    qBC[0] = enorm*nx[0] - etang*nx[1];
    qBC[1] = enorm*nx[1] + etang*nx[0];
    qBC[2] = -ez;

//    qBC[0] = q[0];
//    qBC[1] = q[1];
//    qBC[2] = q[2];

    // B-field
    // RMF
    REAL t = xc[0]; // current time
    REAL Bt = _B0*(1.-exp(-t/_rise))*cos(2.*_pi*_freq*t+_phase);

    qBC[3] = Bt*nx[1];
    qBC[4] =-Bt*nx[0];
    qBC[5] = q[5];

    qBC[6] = -q[6];
    qBC[7] = q[7];
}

// instantiations
template class WxMaxwellRMFAntennaBC<float>;
template class WxMaxwellRMFAntennaBC<double>;
