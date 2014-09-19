#ifndef __wxhyperbolicsrcset__
#define __wxhyperbolicsrcset__

// WarpX hyper includes
#include <wxhyperbolicsrc.h>

// WarpX lib includes
#include <wxobject.h>
#include <wxexcept.h>

template<typename REAL>
class WxHyperbolicSrcSet : public WxObject
{
  public:
/** Create a new object */
    WxHyperbolicSrcSet() {}

/** Dtor */
    virtual ~WxHyperbolicSrcSet();

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
 * Computes source terms from individual sources and assembles them
 * into the full source term for the equation system being solved
 */
    void sourceTerms(REAL *tx, REAL *q, REAL *qaux, REAL *s);

  private:
/** Private to prevent copying */
    WxHyperbolicSrcSet(const WxHyperbolicSrcSet<REAL>&);

/** Private to prevent assignment */
    WxHyperbolicSrcSet<REAL>& operator=(const WxHyperbolicSrcSet<REAL>&);

/** List of sources in system */
    std::vector<WxHyperbolicSrc<REAL>* > _src;
/** Number of equations */
    unsigned _meqn;
};

#endif // __wxhyperbolicsrcset__
