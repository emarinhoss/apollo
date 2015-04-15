// WarpX subsolver includes
#include "wxgridbc.h"

template <typename REAL>
void
WxGridBC<REAL>::setup(const WxCryptSet& wxc, DM dm)
{
  // setup our parent subsolver
  ApSubSolver<REAL>::setup(wxc, dm);

}

template <typename REAL>
WxStepperStatus<REAL>
WxGridBC<REAL>::step(REAL t, REAL dt, Vec in, Vec out)
{
  // applyToArray(dt, in);
  return WxStepperStatus<REAL>();
}

template <typename REAL>
void
WxGridBC<REAL>::applyToArray(REAL *xc, REAL *nx, REAL *q, REAL *qaux, REAL *qBC)
{
    this->applyBC(xc,nx,q,qaux,qBC);
}

// instantiations
template class WxGridBC<float>;
template class WxGridBC<double>;
