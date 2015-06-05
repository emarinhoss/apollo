#ifndef __wxchargesrc__h__
#define __wxchargesrc__h__

// WarpX includes

// WarpX hyperbolic solver includes
#include <wxhyperbolicsrc.h>

// std includes
#include <iostream>
#include <string>

template<class REAL>
class WxChargeSrc : public WxHyperbolicSrc<REAL>
{
public:

    WxChargeSrc()
            : WxHyperbolicSrc<REAL>("charges") {
    }

    void setup(const WxCryptSet& wxc) {
        // call base-class setup function
        WxHyperbolicSrc<REAL>::setup(wxc);

        // read charge and mass of particles
        REAL _q = wxc.template get<REAL>("charge");
        REAL _m = wxc.template get<REAL> ("mass");
        REAL _epsilon0 = wxc.template get<REAL> ("epsilon0");
        REAL _chi = wxc.template get<REAL> ("chi");

        _qbyme = _q*_chi/(_m*_epsilon0);
    }
    
/**
 * Calculate value of source terms at given time/space locations for
 * given input conserved variables.
 *
 * @param n Length of q array
 * @param tx[0] is time and tx[1..3] are spatial location
 * @param q Conserved variables at tx
 * @param qaux Conserved variables at tx
 * @param s Source term
 */
    bool src(unsigned n, REAL *tx, REAL *q, REAL *qaux, REAL *s) {
        // assumes that q is [rho]

        // compute source terms
        s[0] = _qbyme*q[0]; // phi equation

        return true;
    }

private:
    REAL _qbyme; // charge/mass ratio
};

#endif //  __wxchargesrc__h__
