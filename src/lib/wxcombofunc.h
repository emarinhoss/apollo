#ifndef __wxcombofunc__h__
#define __wxcombofunc__h__

// WarpX includes
#include "wxfunction.h"
#include "wxcryptset.h"
#include "wxcreator.h"

// std includes
#include <vector>
#include <string>
#include <iostream>

template<typename INP>
class WxComboFunc : public WxFunction<INP>
{
  public:

/**
 * Size of output array
 */
    unsigned size() const {
      return _nvars;
    }

/**
 * Setup the function
 */
    void setup(const WxCryptSet& wxc);

/**
 * Evaluate the function
 */
    bool func(unsigned n, INP *x, INP *r) {

      unsigned loc = 0;

      // loop over each sub-function....
      typename std::vector<WxFunction<INP>* >::const_iterator i;
      for (i=_funcs.begin(); i!= _funcs.end(); ++i)
      {
        // ... call it ....
        (*i)->func(n, x, r+loc);
        // ... and move the location pointer 
        loc += (*i)->size();
      }

      return true;
    }

  private:
    unsigned _nvars;
    std::vector<WxFunction<INP>* > _funcs;
};

template<typename INP>
void
WxComboFunc<INP>::setup(const WxCryptSet& wxc)
{
  // read list of functions in this combination
  std::vector<WxAny> fnames = wxc.template get<std::vector<WxAny> >("functions");

  _nvars = 0;
  // loop over each function adding it to the list
  typename std::vector<WxAny>::const_iterator i;
  for (i = fnames.begin(); i!= fnames.end(); ++i)
  {
    std::string fname = wx_any_cast<std::string>(*i);
    // find the function's cryptset
    const WxCryptSet& fcs = wxc.getSet( fname );
    // find its kind field
    std::string kind = fcs.template get<std::string>("Kind");        
    // create the function
    WxFunction<INP> *f = WxCreatorMap<WxFunction<INP> >::getNew(kind);
    // setup the function with its cryptset
    f->setup(fcs);
    // accumulate the no of components set by this funtion
    _nvars += f->size();
    // add it to list
    _funcs.push_back(f);
  }
}

#endif //  __wxcombofunc__h__
