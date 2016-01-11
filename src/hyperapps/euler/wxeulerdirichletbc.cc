// WarpX subsolver includes
#include "wxeulerdirichletbc.h"

template <typename REAL>
void
WxEulerDirichletBC<REAL>::setup(const WxCryptSet& wxc, DM dm)
{
  // call base class setup
  WxGridBC<REAL>::setup(wxc, dm);

  _rho = wxc.template get<REAL>("density");
  _u = wxc.template get<REAL>("xvelocity");
  _v = wxc.template get<REAL>("yvelocity");
  _w = wxc.template get<REAL>("zvelocity");
  _gamma = wxc.template get<REAL>("gas_gamma");
  _p = wxc.template get<REAL>("pressure");

}

template <typename REAL>
void
WxEulerDirichletBC<REAL>::applyBC(REAL *xc, REAL *nx, REAL *q, REAL *qaux, REAL *AreaInts, REAL *qBC)
{
    qBC[0] = _rho;
    qBC[1] = _rho*_u;
    qBC[2] = _rho*_v;
    qBC[3] = _rho*_w;
    qBC[4] = _p/(_gamma-1) + 0.5*_rho*(_u*_u+_v*_v+_w*_w);
}

// instantiations
template class WxEulerDirichletBC<float>;
template class WxEulerDirichletBC<double>;
