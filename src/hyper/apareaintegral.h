#ifndef APAREAINTEGRAL_H
#define APAREAINTEGRAL_H

// WarpX includes
#include <wxany.h>
#include <wxobject.h>
#include <wxcryptset.h>

// std includes
#include <string>
#include <vector>

/**
 * Base class to represent area integral terms for use in hyperbolic
 * equations.
 */
template <typename REAL>
class ApAreaIntegral : public WxObject
{
  public:
/**
 * Create a new area integral term with given name
 */
    ApAreaIntegral(const std::string& name);

/**
 *  Setup object using supplied crypset. This method should read in
 *  parameters needed to initialize it.
 *
 * @param wxc Cryptset using which the object is set up.
 */
    virtual void setup(const WxCryptSet& wxc);

/**
 * Calculate value of area integral at given time/space locations
 * for given input conserved variables. This must be provided in
 * derived classes.
 *
 * @param n Length of q array
 * @param tx[0] is time and tx[1..3] are spatial location
 * @param q Conserved variables at tx
 * @param qaux Auxillary variables at tx
 * @param s Source term
 */
    virtual bool integral(unsigned n, REAL *tx, REAL *q, REAL *qaux, REAL *s) = 0;

/**
 * Calculates source term.
 */
    bool compAreaIntegral(REAL *tx, REAL *qfull, REAL *qauxfull, REAL *sfull);

private:
  unsigned _ninp, _nauxinp;
  std::vector<unsigned> _inpIndices, _inpAuxIndices; // input indices
  unsigned _nout;
  std::vector<unsigned> _outIndices; // output indices
  std::vector<REAL> _inpValues, _inpAuxValues, _outValues; // actual input/output data
};

#endif // APAREAINTEGRAL_H
