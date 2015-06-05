#ifndef __wxcurrentsrc__h__
#define __wxcurrentsrc__h__

// WarpX includes

// WarpX hyperbolic solver includes
#include <wxhyperbolicsrc.h>

// std includes
#include <iostream>
#include <string>

template<class REAL>
class WxCurrentSrc : public WxHyperbolicSrc<REAL>
{
public:

    WxCurrentSrc()
            : WxHyperbolicSrc<REAL>("currents") 
      {
      }

    void setup(const WxCryptSet& wxc) {
        // call base-class setup function
        WxHyperbolicSrc<REAL>::setup(wxc);

        // read charge and mass of particles
        _q = wxc.template get<REAL>("charge");
        _m = wxc.template get<REAL> ("mass");
        _epsilon0 = wxc.template get<REAL> ("epsilon0");

        _qbyme = _q/(_m*_epsilon0);
    }
    
/**
 * Calculate value of source terms at given time/space locations for
 * given input conserved variables.
 *
 * @param n Length of q array
 * @param tx[0] is time and tx[1..3] are spatial location
 * @param q Conserved variables at tx
 * @param qaux Auxillary variables
 * @param s Source term
 */
    bool src(unsigned n, REAL *tx, REAL *q, REAL *qaux, REAL *s) {
        // assumes that q is [rho*u, rho*v, rho*w]

        REAL rhou = q[0];
        REAL rhov = q[1];
        REAL rhow = q[2];

        // compute source terms
        s[0] = -_qbyme*rhou;
        s[1] = -_qbyme*rhov;
        s[2] = -_qbyme*rhow;

        return true;
    }

private:
    REAL _q, _m, _epsilon0, _qbyme; // charge, mass and charge/mass ratio
};

#endif //  __wxcurrentsrc__h__
