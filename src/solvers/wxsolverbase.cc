// WarpX solver includes
#include "wxsolverbase.h"

template <typename REAL>
WxSolverBase<REAL>::WxSolverBase(const std::string& name)
  : WxStepper<REAL>(name)
{
}

template <typename REAL>
WxSolverBase<REAL>::~WxSolverBase()
{
}

template <typename REAL>
void
WxSolverBase<REAL>::setRunName(const std::string& runName) 
{
  _runName = runName;
}

template <typename REAL>
std::string
WxSolverBase<REAL>::runName() const 
{
  return _runName;
}

template <typename REAL>
void
WxSolverBase<REAL>::setSolverName(const std::string& name)
{
  _solverName = name;
}

template <typename REAL>
std::string
WxSolverBase<REAL>::getSolverName() const
{
  return _solverName;
}

// instantiations
//template class WxSolverBase<float>;
template class WxSolverBase<double>;
