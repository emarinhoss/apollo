#ifndef APRMFSRC_H
#define APRMFSRC_H

// WarpX includes
#include <wxlogger.h>

// WarpX hyperbolic solver includes
#include <wxhyperbolicsrc.h>

// std includes
#include <iostream>
#include <string>
#include <cmath>

template<class REAL>
class ApRMFSrc : public WxHyperbolicSrc<REAL>
{
  public:

    ApRMFSrc()
      : WxHyperbolicSrc<REAL>("maxwellRMFSrc") {
    }

    void setup(const WxCryptSet& wxc) {
      // call base-class setup function
      WxHyperbolicSrc<REAL>::setup(wxc);

      REAL freq = wxc.template get<REAL>("frequency");
      _B0 = wxc.template get<REAL>("B_rmf");
      _phase = wxc.template get<REAL>("phase");
      _rise = wxc.template get<REAL>("rise_time");
      _r0 = wxc.template get<REAL>("radius");

      _pi = 3.141592653589793;

      _omega = 2*_pi*freq;
    }

/**
 * Calculate value of Braginskii transport terms for the momentum and
 * energy equations at given time/space locations for given input
 * conserved and auxiliary variables.
 *
 * @param n Length of q array
 * @param tx[0] is time and tx[1..3] are spatial location
 * @param q Conserved variables at tx
 * @param qaux Auxiliary variables at tx
 * @param s Source term
 */
    bool src(unsigned n, REAL *tx, REAL *q, REAL *qaux, REAL *s) {

        REAL x = tx[1];
        REAL y = tx[2];
        REAL r = sqrt(x*x+y*y);
        REAL theta = atan(y/x);

        // RMF
        REAL t = tx[0]; // current time
        REAL Bt = _B0*(1.-exp(-t/_rise));

        REAL Br = Bt*sin(_omega*t+_phase);
        REAL Bc = Bt*cos(_omega*t+_phase);

        REAL Ez = r*(_B0*(-exp(-t/_rise))/_rise*cos(_omega*t+_phase)-Bt*_omega*sin(_omega*t+_phase));


        if(r>_r0){
            s[1] = _omega*Br;
            s[2] = _omega*Bc;
            s[0] = _omega*Ez;
        }

      return true;
    }

  private:
    REAL _omega, _B0, _pi, _rise, _phase, _r0;

};

#endif // APRMFSRC_H
