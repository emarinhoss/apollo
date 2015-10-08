#ifndef WXPOINTCURRENT_H
#define WXPOINTCURRENT_H

// WarpX includes

// WarpX hyperbolic solver includes
#include <wxhyperbolicsrc.h>

// std includes
#include <iostream>
#include <string>

template<class REAL>
class WxpointCurrentSrc : public WxHyperbolicSrc<REAL>
{
public:

    WxpointCurrentSrc()
            : WxHyperbolicSrc<REAL>("pointCurrentSrc") {
    }

    void setup(const WxCryptSet& wxc) {
        // call base-class setup function
        WxHyperbolicSrc<REAL>::setup(wxc);

        _xi = wxc.template get<REAL> ("xLocation");
        _yi = wxc.template get<REAL> ("yLocation");
        _r = wxc.template get<REAL> ("radius");
        _j0 = wxc.template get<REAL> ("currentDensity");
        _omega = wxc.template get<REAL> ("frequency");
        _phase = wxc.template get<REAL> ("phase");

        // read charge and mass of particles
        REAL q = wxc.template get<REAL>("charge");
        REAL m = wxc.template get<REAL> ("mass");
        REAL epsilon0 = wxc.template get<REAL> ("epsilon0");

        _qbyme = q/(m*epsilon0);

    }

/**
 * Calculate value of source terms at source at given time/space
 * locations for given input conserved variables.
 *
 * @param n Length of q array
 * @param tx[0] is time and tx[1..3] are spatial location
 * @param q Conserved variables at tx
 * @param qaux Conserved variables at tx
 * @param s Source term
 */
    bool src(unsigned n, REAL *tx, REAL *q, REAL *qaux, REAL *s) {
        // Not dependent on q only. Dependent only on position and time.
        REAL t = tx[0];
        REAL x = tx[1];
        REAL y = tx[2];

        if(sqrt((x-_xi)*(x-_xi)+(y-_yi)*(y-_yi))<_r){
            s[0] = _j0*_qbyme*sin(_omega*t+_phase);
        }
        else{
            s[0] = 0.0;
        }

        return true;
    }

private:
    REAL _xi, _yi, _r, _j0, _qbyme, _omega, _phase; // x/y location of current, radius of "wire", frequency of oscillation, current density
};

#endif // WXPOINTCURRENT_H
