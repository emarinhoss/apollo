#include "apmaxwellrmfantennabc.h"
#include <wxmath.h>

template <typename REAL>
void
WxMaxwellRMFAntennaBC<REAL>::setup(const WxCryptSet& wxc, DM dm)
{
  // call base class setup
  WxGridBC<REAL>::setup(wxc, dm);

  REAL freq = wxc.template get<REAL>("frequency");
  _B0 = wxc.template get<REAL>("B_rmf");
  _phase = wxc.template get<REAL>("phase");
  _rise = wxc.template get<REAL>("rise_time");
  _epsilon_r = wxc.template get<REAL>("relative_permittivity");

  _pi = 3.141592653589793;
  _omega = 2.*_pi*freq;

}

template <typename REAL>
void
WxMaxwellRMFAntennaBC<REAL>::applyBC(REAL *xc, REAL *nx, REAL *q, REAL *qaux, REAL *AreaInts, REAL *qBC)
{
    REAL ex = q[0];
    REAL ey = q[1];

    REAL enorm = ex*nx[0] + ey*nx[1];
    REAL etang = ex*nx[1] - ey*nx[0];
    REAL newEnorm = 1./_epsilon_r*enorm;

    qBC[0] = newEnorm*nx[0] + etang*nx[1];
    qBC[1] = newEnorm*nx[1] - etang*nx[0];

//    qBC[0] = q[0];
//    qBC[1] = q[1];

    // B-field
    // RMF
    REAL t = xc[0]; // current time
    REAL Bt = _B0*(1.-exp(-t/_rise))*cos(_omega*t+_phase);

    qBC[2] = _B0*(-exp(-t/_rise)/_rise*cos(_omega*t+_phase)
                   -(1.-exp(-t/_rise))*sin(_omega*t+_phase)*_omega);
    qBC[3] = Bt*nx[1];
    qBC[4] =-Bt*nx[0];
    qBC[5] = q[5];

    qBC[6] = q[6];
    qBC[7] = q[7];
}

// instantiations
template class WxMaxwellRMFAntennaBC<float>;
template class WxMaxwellRMFAntennaBC<double>;
