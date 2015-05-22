// WarpX subsolver includes
#include "wxgridfixedbc.h"

template <typename REAL>
void
WxGridFixedBC<REAL>::setup(const WxCryptSet& wxc, DM dm)
{
  // call base class setup
  WxGridBC<REAL>::setup(wxc, dm);

  _value = wxc.template get<REAL>("value");

}

template <typename REAL>
void
WxGridFixedBC<REAL>::applyBC(REAL *xc, REAL *nx, REAL *q, REAL *qaux, REAL *qBC)
{
    qBC[0] = _value;
}

// instantiations
template class WxGridFixedBC<float>;
template class WxGridFixedBC<double>;
