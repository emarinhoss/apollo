#include "apsimplifiedrmfbc.h"
#include <wxmath.h>

template <typename REAL>
void
APsimplifiedRMFbc<REAL>::setup(const WxCryptSet& wxc, DM dm)
{
  // call base class setup
  WxGridBC<REAL>::setup(wxc, dm);

  REAL freq = wxc.template get<REAL>("frequency");
  _B0 = wxc.template get<REAL>("B_rmf");
  _rise = wxc.template get<REAL>("rise_time");

  _pi = 3.141592653589793;
  _omega = 2.*_pi*freq;

}

template <typename REAL>
void
APsimplifiedRMFbc<REAL>::applyBC(REAL *xc, REAL *nx, REAL *q, REAL *qaux, REAL *AreaInts, REAL *qBC)
{
    REAL theta = atan(xc[2]/xc[1]);
    REAL r     = sqrt(xc[1]*xc[1]+xc[2]*xc[2]);
    REAL ex = q[0];
    REAL ey = q[1];
    REAL alpha = _pi/6.;

    REAL bz = q[5];
    REAL phi= q[6];
    REAL psi= q[7];

    REAL enorm = ex*nx[0] + ey*nx[1];
    REAL etang =-ex*nx[1] + ey*nx[0];

//    qBC[0] = q[0];
//    qBC[1] = q[1];

    // B-field
    // RMF
    REAL t = xc[0]; // current time
    REAL Bt = 0.5*_B0*(1.-exp(-t/_rise));

    REAL Br =-Bt*(sin(_omega*t+theta)+sin(_omega*t+theta+alpha));
    REAL Bc =-Bt*(cos(_omega*t+theta)+cos(_omega*t+theta+alpha));

    REAL Ez =-r*(0.5*_B0*(-exp(-t/_rise))/_rise*(cos(_omega*t+theta)+cos(_omega*t+theta+alpha))
                                     -Bt*_omega*(sin(_omega*t+theta)+sin(_omega*t+theta+alpha)));

    qBC[0] = enorm*nx[0] + etang*nx[1];
    qBC[1] = enorm*nx[1] - etang*nx[0];
    qBC[2] = Ez;

    qBC[3] = Br;
    qBC[4] = Bc;
    qBC[5] = bz;

    qBC[6] =-phi;
    qBC[7] = psi;
}

// instantiations
template class APsimplifiedRMFbc<float>;
template class APsimplifiedRMFbc<double>;
