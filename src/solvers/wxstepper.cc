// WaprX solver includes
#include <wxstepper.h>

template <typename REAL>
WxStepper<REAL>::WxStepper(const std::string& name)
  : WxObject(name) {
}

template <typename REAL>
WxStepper<REAL>::~WxStepper()
{
}

template <typename REAL>
void
WxStepper<REAL>::setDt(REAL dt)
{
  _dt = dt;
}

template <typename REAL>
void
WxStepper<REAL>::setCurrentTime(REAL tcurr)
{
  _currTime = tcurr;
}

// instantiations
//template class WxStepper<float>;
template class WxStepper<double>;
