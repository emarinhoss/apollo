#include "wxtwofluidrmfbc.h"
#include <wxmath.h>

template <typename REAL>
void
WxTwoFluidRMFBC<REAL>::setup(const WxCryptSet& wxc, DM dm)
{
  // call base class setup
  WxGridBC<REAL>::setup(wxc, dm);

  _omega = wxc.template get<REAL>("frequency");
  _baxial = wxc.template get<REAL>("B_axial");
  _B0 = wxc.template get<REAL>("B_rmf");
  _phase = wxc.template get<REAL>("phase");

}

template <typename REAL>
void
WxTwoFluidRMFBC<REAL>::applyBC(REAL *xc, REAL *nx, REAL *q, REAL *qaux, REAL *qBC)
{
    // electrons
    qBC[0] = q[0];
    qBC[1] = q[1]-2.*(nx[0]*q[1]+nx[1]*q[2])*nx[0];
    qBC[2] = q[2]-2.*(nx[0]*q[1]+nx[1]*q[2])*nx[1];
    qBC[3] = q[3];
    qBC[4] = q[4];

    // ions
    qBC[5] = q[5];
    qBC[6] = q[6]-2.*(nx[0]*q[6]+nx[1]*q[7])*nx[0];
    qBC[7] = q[7]-2.*(nx[0]*q[6]+nx[1]*q[7])*nx[1];
    qBC[8] = q[8];
    qBC[9] = q[9];

    REAL ex = q[10];
    REAL ey = q[11];
    REAL ez = q[12];
    REAL bx = q[13];
    REAL by = q[14];
    REAL bz = q[15];
    REAL phi= q[16];
    REAL psi= q[17];

    // E-field
    REAL enorm = ex*nx[0] + ey*nx[1];
    REAL etang = ex*nx[1] - ey*nx[0];

    qBC[10] = enorm*nx[0] - etang*nx[1];
    qBC[11] = enorm*nx[1] + etang*nx[0];
    qBC[12] = -ez;

//    qBC[10] = ex;
//    qBC[11] = ey;
//    qBC[12] = ez;

    // B-field
    // RMF
    REAL t = xc[0]; // current time
    REAL Bomega_x = _B0*cos(_omega*t+_phase);
    REAL Bomega_y = _B0*sin(_omega*t+_phase);

//    REAL bnorm = -bx*nx[0] - by*nx[1];
//    REAL btang = -bx*nx[1] + by*nx[0];

//    qBC[13] = bnorm*nx[0] - btang*nx[1] + Bomega_x;
//    qBC[14] = bnorm*nx[1] + btang*nx[0] + Bomega_y;
//    qBC[15] = bz;

    qBC[13] = Bomega_x;
    qBC[14] = Bomega_y;
    qBC[15] = _baxial;

    qBC[16] = -phi;
    qBC[17] =  psi;
}

// instantiations
template class WxTwoFluidRMFBC<float>;
template class WxTwoFluidRMFBC<double>;
