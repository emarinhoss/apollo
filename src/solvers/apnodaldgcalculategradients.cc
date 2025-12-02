#include "apnodaldgcalculategradients.h"

template <typename REAL>
void
ApNodalDGcalculateGradients<REAL>::setup(const WxCryptSet& wxc, DM dm)
{
  // setup our parent subsolver
  ApSubSolver<REAL>::setup(wxc, dm);

}

template <typename REAL>
WxStepperStatus<REAL>
ApNodalDGcalculateGradients<REAL>::step(REAL t, REAL dt, Vec in, Vec out)
{
  // applyToArray(dt, in);
  return WxStepperStatus<REAL>();
}

template <typename REAL>
void
ApNodalDGcalculateGradients<REAL>::calculateGradients(wxNodalDGgeometry2D<REAL> *geom, WxCubature2d<REAL> *cub, Vec qk, Vec q_limited)
{
    this->CalcGrad(geom,cub,qk,q_limited);
}

// instantiations
//template class ApNodalDGcalculateGradients<float>;
template class ApNodalDGcalculateGradients<double>;
