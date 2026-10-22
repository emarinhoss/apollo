#include "apemfielsrmfwithharmonics.h"
#include <wxmath.h>

template <typename REAL>
void
ApEMFieldsRMFWithHarmonics<REAL>::setup(const WxCryptSet& wxc, DM dm)
{
  // call base class setup
  WxGridBC<REAL>::setup(wxc, dm);

  REAL freq = wxc.template get<REAL>("frequency");
  _B0 = wxc.template get<REAL>("B_rmf");
  _phase = wxc.template get<REAL>("phase");
  _harmonics = wxc.template get<int>("numberOfHarmonics");
  _a = wxc.template get<REAL>("coilRadius");

  _pi = 3.141592653589793;
  _omega = 2.*_pi*freq;
}

template <typename REAL>
void
ApEMFieldsRMFWithHarmonics<REAL>::applyBC(REAL *xc, REAL *nx, REAL *q, REAL *qaux, REAL *AreaInts, REAL *qBC)
{
    REAL alpha = _pi/4.;    // angular spacing between the two coils at different phases
    REAL t = xc[0];         // current time
    REAL r = sqrt(xc[1]*xc[1]+xc[2]*xc[2]);
    REAL theta = atan(xc[2]/xc[1]);

    REAL Br=0., Bt=0., Ez=0.;

    REAL ex = q[0];
    REAL ey = q[1];

    REAL phi= q[6];
    REAL psi= q[7];

    REAL enorm =  ex*nx[0] + ey*nx[1];
    REAL etang = -ex*nx[1] + ey*nx[0];

    qBC[0] = enorm*nx[0] - etang*nx[1];
    qBC[1] = enorm*nx[1] + etang*nx[0];

    qBC[5] = q[5];

    qBC[6] =-phi;
    qBC[7] = psi;

    // RMF
    for(unsigned j=1; j<_harmonics+1; j++){
        Br += -pow(-1.,j)*cos(0.5*(2.*j-1.)*alpha)/r           *pow(r/_a,2.*j-1.)*cos(_omega*t+pow(-1.,j)*(2.*j-1.)*theta);
        Bt +=             cos(0.5*(2.*j-1.)*alpha)/(2.*j-1.)/_a*pow(r/_a,2.*j-2.)*sin(_omega*t+pow(-1.,j)*(2.*j-1.)*theta);
        Ez +=      _omega*cos(0.5*(2.*j-1.)*alpha)/(2.*j-1.)   *pow(r/_a,2.*j-1.)*cos(_omega*t+pow(-1.,j)*(2.*j-1.)*theta);
    }

    qBC[2] = _B0*_a*Ez;

    qBC[3] = _B0*_a*Bt;
    qBC[4] = _B0*_a*Br;

}

// instantiations
//template class ApEMFieldsRMFWithHarmonics<float>;
template class ApEMFieldsRMFWithHarmonics<double>;
