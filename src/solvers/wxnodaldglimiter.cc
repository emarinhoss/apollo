#include "wxnodaldglimiter.h"

template <typename REAL>
void
WxNodalDGLimiter<REAL>::setup(const WxCryptSet& wxc, DM dm)
{
  // setup our parent subsolver
  ApSubSolver<REAL>::setup(wxc, dm);

}

template <typename REAL>
WxStepperStatus<REAL>
WxNodalDGLimiter<REAL>::step(REAL t, REAL dt, Vec in, Vec out)
{
  // applyToArray(dt, in);
  return WxStepperStatus<REAL>();
}

template <typename REAL>
void
WxNodalDGLimiter<REAL>::applyToVector(wxNodalDGgeometry2D<REAL> *geom, Vec qk, Vec q_limited)
{
    this->applyLimiter(geom,qk,q_limited);
}

// instantiations
template class WxNodalDGLimiter<float>;
template class WxNodalDGLimiter<double>;
