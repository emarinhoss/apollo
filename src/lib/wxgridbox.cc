#include "wxgridbox.h"


template<typename REAL>
WxGridBox<REAL>::WxGridBox(unsigned ndims)
  : WxBox<REAL>(ndims)
{
  for (unsigned i=0; i<ndims; ++i)
    _isDirPeriodic[i] = false;
}

template<typename REAL>
WxGridBox<REAL>::WxGridBox(unsigned dims, REAL *lower, REAL *upper, const WxRange& r, unsigned npad)
  : WxBox<REAL>(dims,lower,upper),  WxGridRange(r, npad)
{
  for (unsigned i=0; i<dims; ++i)
    _isDirPeriodic[i] = false;
}

template<typename REAL>
WxGridBox<REAL>::WxGridBox(const WxGridBox<REAL>& gb)
  : WxBox<REAL>(gb), WxGridRange(gb)
{
  for (unsigned i=0; i<gb.ndims(); ++i)
    _isDirPeriodic[i] = gb._isDirPeriodic[i];
}

template<typename REAL>
WxGridBox<REAL>&
WxGridBox<REAL>::operator=(const WxGridBox<REAL>& gb)
{
  if (this==&gb) return *this;
  // call base class assignment operators
  WxBox<REAL>::operator= (gb);
  WxGridRange::operator= (gb);
  for (unsigned i=0; i<gb.ndims(); ++i)
    _isDirPeriodic[i] = gb._isDirPeriodic[i];
  return *this;
}

template<typename REAL>
void
WxGridBox<REAL>::setup(const WxCryptSet& wxc)
{
  // setup the parent box object
  WxBox<REAL>::setup(wxc);

  unsigned npad = 0;
  if (wxc.has("Pad"))
    npad = wxc.template get<int>("Pad");
  std::vector<WxAny> cells = wxc.template get<std::vector<WxAny> >("Cells");

  unsigned dims = this->ndims();
  int *lower_idx = new int[dims];
  int *upper_idx = new int[dims];
    
  // construct interior range object
  for (unsigned i=0; i < dims; ++i)
  {
    lower_idx[i] = 0;
    upper_idx[i] = wx_any_cast<int>(cells[i]);
  }
  WxRange r = WxRange(dims, lower_idx, upper_idx);

  // construct full range object
  for (unsigned i=0; i < dims; ++i)
  {
    lower_idx[i] = r.lower(i) - npad;
    upper_idx[i] = r.upper(i) + npad;
  }
  WxRange fr = WxRange(dims, lower_idx, upper_idx);

  // set ranges 
  this->setRanges(npad, r, fr);

  // set periodic directions
  for (unsigned i=0; i<dims; ++i)
    _isDirPeriodic[i] = false;

  if (wxc.has("PeriodicDirs")) 
  {
    std::vector<WxAny> periodicDirs = 
      wxc.template get<std::vector<WxAny> >("PeriodicDirs");
    for (unsigned i=0; i<periodicDirs.size(); ++i)
      _isDirPeriodic[wx_any_cast<int>(periodicDirs[i])] = true;
  }

  if (wxc.has("decomp"))
  {
    std::vector<WxAny> decomp = 
      wxc.template get<std::vector<WxAny> >("decomp");
    if (decomp.size() != dims) 
    {
      WxExcept wxe("WxGridBox::setup: Length of 'decomp' array should be ");
      wxe << dims << " but is " << decomp.size() << std::endl;
      throw wxe;
    }
    for (unsigned i=0; i<dims; ++i)
    {
      _decomp.push_back(wx_any_cast<int>(decomp[i]));
    }
    
  }

  delete [] lower_idx;
  delete [] upper_idx;
}        

// instantiations
//template class WxGridBox<float>;
template class WxGridBox<double>;
