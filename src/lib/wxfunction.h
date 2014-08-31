#ifndef __wxfunction__h__
#define __wxfunction__h__

// WarpX includes
#include "wxcryptset.h"
#include "wxobject.h"

// std includes
#include <vector>

template<typename INP, typename OUT=INP>
class WxFunction : public WxObject
{
  public:

    virtual ~WxFunction() {
    }

/**
 * Size of result array. Default is 1
 */
    virtual unsigned size() const {
      return 1; // by default all wxfunctions return scalar
    }

/**
 * Calculate value of function: returns true if function was computed
 * successfully and false otherwise.
 * 
 * @param n Size of input array
 * @param x Input array of size 'n' of type INP
 * @param r Result array of size this->size() of type OUT
 */
    virtual bool func(unsigned n, INP *x, OUT *r) = 0;
};

#endif // __wxfunction__h__
