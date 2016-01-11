#include "wxphmaxwellconductingwallbc.h"

template <typename REAL>
void
WxPHMaxwellConductingWallBC<REAL>::setup(const WxCryptSet& wxc, DM dm)
{
  // call base class setup
  WxGridBC<REAL>::setup(wxc, dm);

}

template <typename REAL>
void
WxPHMaxwellConductingWallBC<REAL>::applyBC(REAL *xc, REAL *nx, REAL *q, REAL *qaux, REAL *AreaInts, REAL *qBC)
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

    REAL bnorm = -bx*nx[0] - by*nx[1];
    REAL btang = -bx*nx[1] + by*nx[0];

    qBC[3] = bnorm*nx[0] - btang*nx[1];
    qBC[4] = bnorm*nx[1] + btang*nx[0];
    qBC[5] = bz;

    qBC[6] = -phi;
    qBC[7] =  psi;
}

// instantiations
template class WxPHMaxwellConductingWallBC<float>;
template class WxPHMaxwellConductingWallBC<double>;
