#ifndef WXTRANSVERSEMAGNETICSRC_H
#define WXTRANSVERSEMAGNETICSRC_H

// WarpX hyperbolic solver includes
#include <wxhyperbolicsrc.h>

// std includes
#include <iostream>
#include <string>
#include <cmath>

template<class REAL>
class WxTransverseMagneticSrc : public WxHyperbolicSrc<REAL>
{
  public:

    WxTransverseMagneticSrc()
            : WxHyperbolicSrc<REAL>("tranverseMagneticSrc") {
    }

    void setup(const WxCryptSet& wxc) {
      // call base-class setup function
      WxHyperbolicSrc<REAL>::setup(wxc);

    }

/**
 * Calculate the radial source terms at given
 * time/space locations for given input conserved variables.
 *
 * @param n Length of q array
 * @param tx[0] is time and tx[1..3] are spatial location
 * @param q Conserved variables at tx
 * @param qaux Conserved variables at tx
 * @param s Source term
 */
    bool src(unsigned n, REAL *tx, REAL *q, REAL *qaux, REAL *s) {
      // assumes that q is [Bz]

      REAL t  = tx[0], x = tx[1], y = tx[2];

      s[0] = exp(t)*(x*(1-x)*y*(1-y)+2.*x*(1-x)+2*y*(1-y));

      return true;
    }

  private:

};


#endif // WXTRANSVERSEMAGNETICSRC_H
