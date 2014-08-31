#ifndef __wxexprfunc__h__
#define __wxexprfunc__h__

// WarpX includes
#include <wxfunction.h>
#include <wxcryptset.h>
#include <wxany.h>
#include <wxexcept.h>
#include <wxexprparser.h>

// std includes
#include <iostream>
#include <string>
#include <vector>
#include <sstream>

/**
 * Sets the output to a specified expr
 */
class WxExprFunc : public WxFunction<double>
{
  public:
/**
 * Size of output array
 *
 * @return size of output array
 */
    unsigned size() const {
      return _nres;
    }

/** Setup the function */
    void setup(const WxCryptSet& wxc);

/** Evaluate function */
    bool func(unsigned n, double *tx, double *r);

  private:
/** No of result variables added */
    unsigned _nres;
/** Expression parser */
    WxExprParser _parser;
/** Result codes of expressions added */
    std::vector<unsigned> _resCodes;
};

#endif //__wxexprfunc__h__
