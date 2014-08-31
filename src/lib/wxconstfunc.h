#ifndef __wxconstfunc__h__
#define __wxconstfunc__h__

// WarpX includes
#include "wxfunction.h"
#include "wxcryptset.h"

// std includes
#include <iostream>
#include <string>

/**
 * Sets the output to a specified constant
 */
template<typename INP, typename OUT=INP>
class WxConstFunc : public WxFunction<INP, OUT>
{
  public:

/**
 * Size of output array
 */
    unsigned size() const {
        return _nvars;
    }

    void setup(const WxCryptSet& wxc) {
        // no of variables to set
        if (wxc.has("nvars"))
            _nvars = wxc.template get<int>("nvars");
        else
            _nvars = 1;

        if (wxc.has("value"))
            _value = wxc.template get<OUT>("value");
        else
            _value = 0;
    }

    bool func(unsigned n, INP *x, OUT *r) {
        for (unsigned i=0; i<_nvars; ++i)
            r[i] = _value;
        return true;
    }

  private:
    unsigned _nvars;
    OUT _value;
};

#endif //  __wxconstfunc__h__
