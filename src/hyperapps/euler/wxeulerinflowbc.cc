// WarpX subsolver includes
#include "wxeulerinflowbc.h"
#include <wxmath.h>

template <typename REAL>
void
WxEulerInflowBC<REAL>::setup(const WxCryptSet& wxc, DM dm)
{
  // call base class setup
  WxGridBC<REAL>::setup(wxc, dm);

  _Mach = wxc.template get<REAL>("Mach_number");
  _angle = wxc.template get<REAL>("Flow_angle");
  _gamma = wxc.template get<REAL>("gas_gamma");
  _pres = wxc.template get<REAL>("Pressure");
  _rho = wxc.template get<REAL>("density");
}

template <typename REAL>
void
WxEulerInflowBC<REAL>::applyBC(REAL *xc, REAL *nx, REAL *q, REAL *qaux, REAL *AreaInts, REAL *qBC)
{
    REAL pi = 3.141592653589793;

    REAL sound = sqrt(_gamma*_pres/_rho);
    REAL u = _Mach*sound*cos(_angle*pi/180.);
    REAL v = _Mach*sound*sin(_angle*pi/180.);

    qBC[0] = _rho;
    qBC[1] = _rho*u;
    qBC[2] = _rho*v;
    qBC[3] = 0.0;
    qBC[4] = _pres/(_gamma-1) + 0.5*_rho*(u*u+v*v);
}

// instantiations
template class WxEulerInflowBC<float>;
template class WxEulerInflowBC<double>;
