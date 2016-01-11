#ifndef APAREAINTEGRALSET_H
#define APAREAINTEGRALSET_H

// WarpX hyper includes
#include "apareaintegral.h"

// WarpX lib includes
#include <wxobject.h>
#include <wxexcept.h>

template<typename REAL>
class ApAreaIntegralSet : public WxObject
{
  public:
/** Create a new object */
    ApAreaIntegralSet() {}

/** Dtor */
    virtual ~ApAreaIntegralSet();

/**
 * Set number of equations
 *
 * @param meqn Number of equations
 */
    void setNumEqns(unsigned meqn);

/**
 * Setup object using supplied crypset. This method should read in
 * parameters needed to initialize it.
 *
 * @param wxc Cryptset using which the object is set up.
 */
    virtual void setup(const WxCryptSet& wxc);

/**
 * Computes area integrals from individual integrals and assembles them
 * into the full integral term for the equation system being solved
 */
    void areaTerms(REAL *tx, REAL *q, REAL *qaux, REAL *s);

  private:
/** Private to prevent copying */
    ApAreaIntegralSet(const ApAreaIntegralSet<REAL>&);

/** Private to prevent assignment */
    ApAreaIntegralSet<REAL>& operator=(const ApAreaIntegralSet<REAL>&);

/** List of sources in system */
    std::vector<ApAreaIntegral<REAL>* > _src;
/** Number of equations */
    unsigned _meqn;
};

#endif // APAREAINTEGRALSET_H
