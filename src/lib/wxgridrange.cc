#include "wxgridrange.h"

WxGridRange::WxGridRange(const WxRange& r, unsigned npad)
  : _npad(npad), _r(r)
{
  unsigned rank = r.ndims();
  int *lower = new int[rank];
  int *upper = new int[rank];

  for (unsigned i=0; i<rank; ++i)
  {
    lower[i] = r.lower(i) - _npad;
    upper[i] = r.upper(i) + _npad;
  }
  _fr = WxRange(rank, lower, upper);
    
  delete [] lower;
  delete [] upper;
}

WxGridRange::WxGridRange(const WxGridRange& wr)
{
  _npad = wr._npad;
  _r = wr._r;
  _fr = wr._fr;
}

WxGridRange&
WxGridRange::operator=(const WxGridRange& wr)
{
  if (this==&wr) return *this;
  _npad = wr._npad;
  _r = wr._r;
  _fr = wr._fr;
  return *this;
}

void
WxGridRange::setRanges(unsigned npad, const WxRange& r, const WxRange& fr)
{
  _npad = npad;
  _r = r;
  _fr = fr;
}

WxRange 
WxGridRange::lowerShell(unsigned dim) const 
{
  return lowerLayer(dim, _npad);
}

WxRange 
WxGridRange::upperShell(unsigned dim) const
{
  return upperLayer(dim, _npad);
}

WxRange 
WxGridRange::lowerLayer(unsigned dim, unsigned nl) const 
{
  unsigned n = this->ndims();
  int lower[16];
  int upper[16];
        
  // set range to that of grid interior
  for (unsigned i=0; i<n; ++i)
  {
    lower[i] = _r.lower(i);
    upper[i] = _r.upper(i);
  }
  // adjust the 'dim' ranges
  upper[dim] = lower[dim] + nl;
  WxRange r(n, lower, upper);

  return r;
}

WxRange 
WxGridRange::upperLayer(unsigned dim, unsigned nl) const
{
  unsigned n = this->ndims();
  int lower[16];
  int upper[16];
        
  // set range to that of grid interior
  for (unsigned i=0; i<n; ++i)
  {
    lower[i] = _r.lower(i);
    upper[i] = _r.upper(i);
  }
  // adjust the 'dim' ranges
  lower[dim] = upper[dim] - nl;
  WxRange r(n, lower, upper);

  return r;
}

WxRange 
WxGridRange::lowerLayerWithGhost(unsigned dim, unsigned nl) const 
{
  unsigned n = this->ndims();
  int lower[16];
  int upper[16];
        
  // set range to that of full grid
  for (unsigned i=0; i<n; ++i)
  {
    lower[i] = _fr.lower(i);
    upper[i] = _fr.upper(i);
  }
  // adjust the 'dim' ranges
  lower[dim] = _r.lower(dim);
  upper[dim] = lower[dim] + nl;
  WxRange r(n, lower, upper);

  return r;
}

WxRange 
WxGridRange::upperLayerWithGhost(unsigned dim, unsigned nl) const
{
  unsigned n = this->ndims();
  int lower[16];
  int upper[16];
        
  // set range to that of grid interior
  for (unsigned i=0; i<n; ++i)
  {
    lower[i] = _fr.lower(i);
    upper[i] = _fr.upper(i);
  }
  // adjust the 'dim' ranges
  upper[dim] = _r.upper(dim);
  lower[dim] = upper[dim] - nl;
  WxRange r(n, lower, upper);

  return r;
}

WxRange 
WxGridRange::lowerGhost(unsigned dim) const 
{
  unsigned n = this->ndims();
  int lower[16];
  int upper[16];
        
  // set lower
  for (unsigned i=0; i<n; ++i)
  {
    lower[i] = _fr.lower(i);
    upper[i] = _r.lower(i);
  }
  WxRange r(n, lower, upper);

  return r;
}

WxRange 
WxGridRange::upperGhost(unsigned dim) const
{
  unsigned n = this->ndims();
  int lower[16];
  int upper[16];
        
  // set lower
  for (unsigned i=0; i<n; ++i)
  {
    lower[i] = _r.upper(i);
    upper[i] = _fr.upper(i);
  }
  WxRange r(n, lower, upper);

  return r;
}
