// WarpX subsolver includes
#include "wxgridfixedbc.h"

template <typename REAL>
void
WxGridFixedBC<REAL>::setup(const WxCryptSet& wxc, DM dm)
{
  // call base class setup
  WxGridBC<REAL>::setup(wxc, dm);

  _value = wxc.template get<REAL>("value");

  // get the label where the BC is to be applied
  _label = wxc.get<std::string>("label");

}

template <typename REAL>
void
WxGridFixedBC<REAL>::applyBC(WxpDGGeometry<REAL> quad, REAL dt, Vec inOut)
{

}

// instantiations
template class WxGridFixedBC<float>;
template class WxGridFixedBC<double>;
