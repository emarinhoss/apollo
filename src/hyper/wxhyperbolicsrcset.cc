// WarpX hyper includes
#include <wxhyperbolicsrcset.h>

// WarpX lib includes
#include <wxlogger.h>
#include <wxlogstream.h>

template<typename REAL>
void
WxHyperbolicSrcSet<REAL>::setNumEqns(unsigned meqn)
{
  _meqn = meqn;
}

template<typename REAL>
void
WxHyperbolicSrcSet<REAL>::setup(const WxCryptSet& wxc)
{
  // get hold of stream to log debug messages
  WxLogStream infoStrm
    = WxLogger::get("apollo-root.console")->getInfoStream();

  // source terms to use, if any
  if (wxc.has("Sources"))
  {
    std::vector<WxAny> src = wxc.template get<std::vector<WxAny> >("Sources");
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
      WxHyperbolicSrc<REAL> *s = WxCreatorMap<WxHyperbolicSrc<REAL> >::getNew(kind);
      // set it up and add it to set
      s->setup(scs);
      _src.push_back(s);
            
      infoStrm << "Source " << srcName << " is of kind " << kind << std::endl;
    }
  }
}

template <typename REAL>
WxHyperbolicSrcSet<REAL>::~WxHyperbolicSrcSet()
{
  WxHyperbolicSrc<REAL>* tempsrc;
  for (int i=0; i<(int)_src.size(); ++i)
  {
    tempsrc = _src[i];
    delete tempsrc;
  }
}

template<typename REAL>
void
WxHyperbolicSrcSet<REAL>::sourceTerms(REAL *tx, REAL *q, REAL *qaux, REAL *s)
{
  // initialize source to 0.0
  for (unsigned m=0; m<_meqn; ++m)
    s[m] = 0.0;

  // loop over each source, adding its contribution to the total source
  typename std::vector<WxHyperbolicSrc<REAL>* >::const_iterator i;
  for (i=_src.begin(); i!=_src.end(); ++i)
    (*i)->compSource(tx, q, qaux, s);
}

// instantiations
template class WxHyperbolicSrcSet<float>;
template class WxHyperbolicSrcSet<double>;
