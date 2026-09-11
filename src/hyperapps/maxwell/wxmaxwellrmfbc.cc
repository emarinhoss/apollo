#include "wxmaxwellrmfbc.h"

#include <wxmath.h>

template <typename REAL>
void
WxMaxwellRMFBC<REAL>::setup(const WxCryptSet& wxc, DM dm)
{
  // call base class setup
  WxGridBC<REAL>::setup(wxc, dm);

  _omega = wxc.template get<REAL>("frequency");
  _B0 = wxc.template get<REAL>("B_rmf");
  _phase = wxc.template get<REAL>("phase");

}

template <typename REAL>
void
WxMaxwellRMFBC<REAL>::applyBC(REAL *xc, REAL *nx, REAL *q, REAL *qaux, REAL *AreaInts, REAL *qBC)
{

    REAL ex = q[0];
    REAL ey = q[1];
    REAL ez = q[2];
    REAL bx = q[3];
    REAL by = q[4];
    REAL bz = q[5];
    REAL phi= q[6];
    REAL psi= q[7];

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
    REAL Bomega_x = _B0*cos(_omega*t+_phase);
    REAL Bomega_y = _B0*sin(_omega*t+_phase);

    qBC[3] = Bomega_x;
    qBC[4] = Bomega_y;
    qBC[5] = q[5];

    qBC[6] = -q[6];
    qBC[7] = q[7];
}

// instantiations
//template class WxMaxwellRMFBC<float>;
template class WxMaxwellRMFBC<double>;
