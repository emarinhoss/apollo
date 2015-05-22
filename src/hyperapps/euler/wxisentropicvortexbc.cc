// WarpX subsolver includes
#include "wxisentropicvortexbc.h"

template <typename REAL>
void
WxIsentropicVortexBC<REAL>::setup(const WxCryptSet& wxc, DM dm)
{
  // call base class setup
  WxGridBC<REAL>::setup(wxc, dm);

  _vo = wxc.template get<REAL>("v0");
  _uo = wxc.template get<REAL>("u0");
  _gamma = wxc.template get<REAL>("gas_gamma");
  _beta = wxc.template get<REAL>("beta");
  _xo = wxc.template get<REAL>("x0");
  _yo = wxc.template get<REAL>("y0");

}

template <typename REAL>
void
WxIsentropicVortexBC<REAL>::applyBC(REAL *xc, REAL *nx, REAL *q, REAL *qaux, REAL *qBC)
{
    REAL t = xc[0];
    REAL x = xc[1];
    REAL y = xc[2];
    REAL pi = 3.141592653589793;

    REAL xmut = x-_uo*t;
    REAL ymvt = y-_vo*t;
    REAL r = sqrt((xmut-_xo)*(xmut-_xo)+(ymvt-_yo)*(ymvt-_yo));
    REAL u = _uo - _beta*exp(1.-r*r)*(ymvt-_yo)/(2.*pi);
    REAL v = _vo + _beta*exp(1.-r*r)*(xmut-_xo)/(2.*pi);
    REAL base = 1. - ( (_gamma-1.)*_beta*_beta*exp(2.*(1.-r*r))/(16.*_gamma*pi*pi) );
    REAL rho1 = pow(base,1./(_gamma-1));
    REAL p1 = pow(rho1,_gamma);

    qBC[0] = rho1;
    qBC[1] = rho1*u;
    qBC[2] = rho1*v;
    qBC[3] = 0.0;
    qBC[4] = p1/(_gamma-1) + 0.5*rho1*(u*u+v*v);
}

// instantiations
template class WxIsentropicVortexBC<float>;
template class WxIsentropicVortexBC<double>;
