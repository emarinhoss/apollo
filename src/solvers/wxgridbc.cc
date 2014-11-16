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
WxGridBC<REAL>::step(REAL dt, Vec in, Vec out)
{
  // applyToArray(dt, in);
  return WxStepperStatus<REAL>();
}

template <typename REAL>
void
WxGridBC<REAL>::applyToArray(WxpDGGeometry<REAL> quad, REAL dt, Vec inOut)
{
    this->applyBC(quad, dt, inOut);
}

// instantiations
template class WxGridBC<float>;
template class WxGridBC<double>;
