#include "apphmaxwelldielectricbc.h"

template <typename REAL>
void
WxPHMaxwellDielectricBC<REAL>::setup(const WxCryptSet& wxc, DM dm)
{
  // call base class setup
  WxGridBC<REAL>::setup(wxc, dm);

  _epsilon_r = wxc.template get<REAL>("relative_permittivity");
  _mu_r = wxc.template get<REAL>("relative_permeability");

}

template <typename REAL>
void
WxPHMaxwellDielectricBC<REAL>::applyBC(REAL *xc, REAL *nx, REAL *q, REAL *qaux, REAL *AreaInts, REAL *qBC)
{
    REAL ex = q[0];
    REAL ey = q[1];
    REAL ez = q[2];
    REAL bx = q[3];
    REAL by = q[4];
    REAL bz = q[5];
    REAL phi= q[6];
    REAL psi= q[7];

    // Electric Field BCs
    REAL enorm = ex*nx[0] + ey*nx[1];
    REAL etang = ex*nx[1] - ey*nx[0];

    REAL newEnorm = 1./_epsilon_r*enorm;

    qBC[0] = newEnorm*nx[0] + etang*nx[1];
    qBC[1] = newEnorm*nx[1] - etang*nx[0];
    qBC[2] = ez;

    // Magnetic Field BCs
    REAL bnorm = bx*nx[0] + by*nx[1];
    REAL btang = bx*nx[1] - by*nx[0];

    REAL newBtang = _mu_r*btang;

    qBC[3] = bnorm*nx[0] + newBtang*nx[1];
    qBC[4] = bnorm*nx[1] - newBtang*nx[0];
    qBC[5] = _mu_r*bz;

    // Corrections
    qBC[6] =-phi;
    qBC[7] = psi;

}

// instantiations
//template class WxPHMaxwellDielectricBC<float>;
template class WxPHMaxwellDielectricBC<double>;
