// WarpX hyper includes
#include "apareaintegralset.h"

// WarpX lib includes
#include <wxlogger.h>
#include <wxlogstream.h>

template<typename REAL>
void
ApAreaIntegralSet<REAL>::setNumEqns(unsigned meqn)
{
  _meqn = meqn;
}

template<typename REAL>
void
ApAreaIntegralSet<REAL>::setup(const WxCryptSet& wxc)
{
  // get hold of stream to log debug messages
  WxLogStream infoStrm
    = WxLogger::get("apollo-root.console")->getInfoStream();

  // source terms to use, if any
  if (wxc.has("AreaIntegrals"))
  {
    std::vector<WxAny> src = wxc.template get<std::vector<WxAny> >("AreaIntegrals");
    // loop over each source adding it to list of sources being solved
    std::vector<WxAny>::const_iterator i;
    for (i=src.begin(); i!=src.end(); ++i)
    {
      // name of cryptset
      std::string srcName = wx_any_cast<std::string>(*i);
      // find the cryptset for this source
      const WxCryptSet scs = wxc.getSet(srcName);

      // find its Kind and create new source object
      std::string kind = scs.template get<std::string>("Kind");
      ApAreaIntegral<REAL> *s = WxCreatorMap<ApAreaIntegral<REAL> >::getNew(kind);
      // set it up and add it to set
      s->setup(scs);
      _src.push_back(s);

      infoStrm << "Area Integral " << srcName << " is of kind " << kind << std::endl;
    }
  }
}

template <typename REAL>
ApAreaIntegralSet<REAL>::~ApAreaIntegralSet()
{
  ApAreaIntegral<REAL>* tempsrc;
  for (int i=0; i<(int)_src.size(); ++i)
  {
    tempsrc = _src[i];
    delete tempsrc;
  }
}

template<typename REAL>
void
ApAreaIntegralSet<REAL>::areaTerms(REAL *tx, REAL *q, REAL *qaux, REAL *s)
{
  // initialize source to 0.0
  for (unsigned m=0; m<_meqn; ++m)
    s[m] = 0.0;

  // loop over each source, adding its contribution to the total source
  typename std::vector<ApAreaIntegral<REAL>* >::const_iterator i;
  for (i=_src.begin(); i!=_src.end(); ++i)
    (*i)->compAreaIntegral(tx, q, qaux, s);
}

// instantiations
//template class ApAreaIntegralSet<float>;
template class ApAreaIntegralSet<double>;
