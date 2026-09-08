#include "aptwofluidrmfantennabc.h"
#include <wxmath.h>

template <typename REAL>
void
WxTwoFluidRMFAntennaBC<REAL>::setup(const WxCryptSet& wxc, DM dm)
{
  // call base class setup
  WxGridBC<REAL>::setup(wxc, dm);

  REAL freq = wxc.template get<REAL>("frequency");
  // _baxial is used by the flux-conserver expression in applyBC() and was
  // never read here, so that expression ran on an uninitialised member.
  _baxial = wxc.template get<REAL>("B_axial");
  _B0 = wxc.template get<REAL>("B_rmf");
  _phase = wxc.template get<REAL>("phase");
  _rise = wxc.template get<REAL>("rise_time");
  _epsilon_r = wxc.template get<REAL>("relative_permittivity");
  _mu_r = wxc.template get<REAL>("relative_permeability");
  _a = wxc.template get<REAL>("plasma_radius");
  _b = wxc.template get<REAL>("flux_conserver_radius");

  _pi = 3.141592653589793;
  _omega = 2.*_pi*freq;

}

template <typename REAL>
void
WxTwoFluidRMFAntennaBC<REAL>::applyBC(REAL *xc, REAL *nx, REAL *q, REAL *qaux, REAL *AreaInts, REAL *qBC)
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

    REAL enorm = ex*nx[0] + ey*nx[1];
    REAL etang = ex*nx[1] - ey*nx[0];
    REAL newEnorm = 1./_epsilon_r*enorm;

    // RMF
    REAL t = xc[0]; // current time
    REAL Bt = _B0*(1.-exp(-t/_rise))*cos(_omega*t+_phase);

    REAL Ez_rmf = _B0*(-exp(-t/_rise)/_rise*cos(_omega*t+_phase)
                       -(1.-exp(-t/_rise))*sin(_omega*t+_phase)*_omega);

    // Area integral \int B_z \cdot dA
    REAL intBzda = AreaInts[15];
    REAL newBz = _b*_b*_baxial/(_b*_b-_a*_a)-intBzda/(_b*_b-_a*_a)/_pi;

    REAL bnorm = bx*nx[0] + by*nx[1];
    REAL btang = bx*nx[1] - by*nx[0];

    REAL newBtang = _mu_r*btang;

    qBC[10] = newEnorm*nx[0] + etang*nx[1];
    qBC[11] = newEnorm*nx[1] - etang*nx[0];
    qBC[12] = Ez_rmf;

    qBC[13] = Bt*nx[1];
    qBC[14] =-Bt*nx[0];
    qBC[15] = newBz;

    qBC[16] = phi;
    qBC[17] = psi;
}

// instantiations
//template class WxTwoFluidRMFAntennaBC<float>;
template class WxTwoFluidRMFAntennaBC<double>;
