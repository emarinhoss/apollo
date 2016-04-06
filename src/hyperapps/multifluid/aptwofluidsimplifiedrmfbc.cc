#include "aptwofluidsimplifiedrmfbc.h"
#include <wxmath.h>

template <typename REAL>
void
APTwoFluidSimplifiedRMFBC<REAL>::setup(const WxCryptSet& wxc, DM dm)
{
  // call base class setup
  WxGridBC<REAL>::setup(wxc, dm);

  REAL freq = wxc.template get<REAL>("frequency");
  _baxial = wxc.template get<REAL>("B_axial");
  _B0 = wxc.template get<REAL>("B_rmf");
  _phase = wxc.template get<REAL>("phase");
  _rise = wxc.template get<REAL>("rise_time");
  _a = wxc.template get<REAL>("plasma_radius");
  _b = wxc.template get<REAL>("flux_conserver_radius");

  _pi = 3.141592653589793;

  _omega = 2*_pi*freq;


}

template <typename REAL>
void
APTwoFluidSimplifiedRMFBC<REAL>::applyBC(REAL *xc, REAL *nx, REAL *q, REAL *qaux, REAL *AreaInts, REAL *qBC)
{
    REAL x = xc[1];
    REAL y = xc[2];
    REAL r = sqrt(x*x+y*y);
    REAL theta = atan(y/x);

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

    REAL phi= q[16];
    REAL psi= q[17];

    // Area integral \int B_z \cdot dA
    REAL intBzda = AreaInts[15];

    // E-field
    REAL enorm = ex*nx[0] + ey*nx[1];
    REAL etang = ex*nx[1] - ey*nx[0];

    qBC[10] = enorm*nx[0] - etang*nx[1];
    qBC[11] = enorm*nx[1] + etang*nx[0];

//    qBC[10] = ex;
//    qBC[11] = ey;

    // B-field
    // RMF
    REAL t = xc[0]; // current time
    REAL Bt = 0.5*_B0*(1.-exp(-t/_rise));

    REAL Br =Bt*sin(_omega*t);
    REAL Bc =Bt*cos(_omega*t);

    REAL Ez =r*(0.5*_B0*(-exp(-t/_rise))/_rise*cos(_omega*t)-Bt*_omega*sin(_omega*t));

    qBC[12] = Ez;
    qBC[13] = Br;
    qBC[14] = Bc;

//    REAL testvalue = 1.413716694115407e-05;
    REAL newBz = _b*_b*_baxial/(_b*_b-_a*_a)-intBzda/(_b*_b-_a*_a)/_pi;
    qBC[15] = newBz;
    qBC[16] = -phi;
    qBC[17] =  psi;
}

// instantiations
template class APTwoFluidSimplifiedRMFBC<float>;
template class APTwoFluidSimplifiedRMFBC<double>;
