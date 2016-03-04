#include "aptwofluiddielectricbc.h"

template <typename REAL>
void
WxTwoFluidDielectricBC<REAL>::setup(const WxCryptSet& wxc, DM dm)
{
  // call base class setup
  WxGridBC<REAL>::setup(wxc, dm);

  _epsilon_r = wxc.template get<REAL>("relative_permittivity");
  _mu_r = wxc.template get<REAL>("relative_permeability");

}

template <typename REAL>
void
WxTwoFluidDielectricBC<REAL>::applyBC(REAL *xc, REAL *nx, REAL *q, REAL *qaux, REAL *AreaInts, REAL *qBC)
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

    // Electric Field BCs
    REAL enorm = ex*nx[0] + ey*nx[1];
    REAL etang = ex*nx[1] - ey*nx[0];

    REAL newEnorm = 1./_epsilon_r*enorm;

    qBC[10] = newEnorm*nx[0] + etang*nx[1];
    qBC[11] = newEnorm*nx[1] - etang*nx[0];
    qBC[12] = ez;

    // Magnetic Field BCs
    REAL bnorm = bx*nx[0] + by*nx[1];
    REAL btang = bx*nx[1] - by*nx[0];

    REAL newBtang = _mu_r*btang;

    qBC[13] = bnorm*nx[0] + newBtang*nx[1];
    qBC[14] = bnorm*nx[1] - newBtang*nx[0];
    qBC[15] = _mu_r*bz;

    // Corrections
    qBC[16] = phi;
    qBC[17] = psi;

}

// instantiations
template class WxTwoFluidDielectricBC<float>;
template class WxTwoFluidDielectricBC<double>;
