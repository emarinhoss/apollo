// WarpX hyper includes
#include <wxhyperbolicsrc.h>

template<typename REAL>
WxHyperbolicSrc<REAL>::WxHyperbolicSrc(const std::string& name)
  : WxObject(name) {
}

template<typename REAL>
void
WxHyperbolicSrc<REAL>::setup(const WxCryptSet& wxc)
{
  _ninp = 0;
  if (wxc.has("InpRange"))
  {
    // read input range
    std::vector<WxAny> inp = 
      wxc.template get<std::vector<WxAny> >("InpRange");
    
    _ninp = inp.size(); // these many components of input array will be used
    _inpIndices.resize(_ninp); 
    _inpValues.resize(_ninp);

    // loop, setting input indices to use
    unsigned j=0;
    typename std::vector<WxAny>::const_iterator i;
    for (i=inp.begin(); i!=inp.end(); ++i)
      _inpIndices[j++] = wx_any_cast<int>(*i);
  }

  _nauxinp = 0;
  if (wxc.has("AuxInpRange"))
  {
    // read aux input range
    std::vector<WxAny> inp = 
      wxc.template get<std::vector<WxAny> >("AuxInpRange");
    _nauxinp = inp.size(); // these many components of input array will be used
    _inpAuxIndices.resize(_nauxinp); 
    _inpAuxValues.resize(_nauxinp);

    // loop, setting input indices to use
    unsigned j=0;
    typename std::vector<WxAny>::const_iterator i;
    for (i=inp.begin(); i!=inp.end(); ++i)
      _inpAuxIndices[j++] = wx_any_cast<int>(*i);
  }

  // read output range
  std::vector<WxAny> out =
    wxc.template get<std::vector<WxAny> >("OutRange");

  _nout = out.size(); // these many components of output array will be used
  _outIndices.resize(_nout); 
  _outValues.resize(_nout);

  // loop, setting output indices to set
  unsigned j=0;
  typename std::vector<WxAny>::const_iterator i;
  for (i=out.begin(); i!=out.end(); ++i)
    _outIndices[j++] = wx_any_cast<int>(*i);
}

template<typename REAL>
bool
WxHyperbolicSrc<REAL>::compSource(REAL *tx, REAL *qfull, REAL *qauxfull, REAL *sfull) 
{
  // copy appropriate components into input array
  for (unsigned i=0; i<_ninp; ++i)
    _inpValues[i] = qfull[_inpIndices[i]];
  
   for (unsigned i=0; i<_nauxinp; ++i)
     _inpAuxValues[i] = qauxfull[_inpAuxIndices[i]];
  
  // make call to compute source term
  bool res = this->src(_ninp, tx, &_inpValues[0], &_inpAuxValues[0], &_outValues[0]);
  
  // add appropriate components from output array
  for (unsigned i=0; i<_nout; ++i)
    sfull[_outIndices[i]] += _outValues[i];
  
  return res;
}

// instantiations
//template class WxHyperbolicSrc<float>;
template class WxHyperbolicSrc<double>;
