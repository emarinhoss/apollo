#ifndef __wxconstvecfunc__h__
#define __wxconstvecfunc__h__

// WarpX includes
#include "wxfunction.h"
#include "wxcryptset.h"
#include "wxany.h"

// std includes
#include <iostream>
#include <string>
#include <vector>

/**
 * Sets the output to a specified const
 */
template<typename INP, typename OUT=INP>
class WxConstVecFunc : public WxFunction<INP, OUT>
{
  public:

/**
 * Size of output array
 */
    unsigned size() const {
        return _nvars;
    }

    void setup(const WxCryptSet& wxc) {
        // read list of values to set
        std::vector<WxAny> vals = 
            wxc.template get<std::vector<WxAny> >("values");
        _nvars = vals.size(); // these many values are to be set
        _value = new INP[_nvars];
        // loop, setting value list
        unsigned j=0;
        typename std::vector<WxAny>::const_iterator i;
        for (i=vals.begin(); i!= vals.end(); ++i)
            _value[j++] = wx_any_cast<INP>(*i);
    }

    bool func(unsigned n, INP *x, OUT *r) {
        for (unsigned i=0; i<_nvars; ++i)
            r[i] = _value[i];
        return true;
    }

  private:
    unsigned _nvars;
    INP *_value;
};

#endif // __wxconstvecfunc__h__
