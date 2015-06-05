#ifndef WXPULSETRAINSRC_H
#define WXPULSETRAINSRC_H

// WarpX hyperbolic solver includes
#include <wxhyperbolicsrc.h>

// std includes
#include <iostream>
#include <string>
#include <cmath>

template<class REAL>
class WxPulseTrainSrc : public WxHyperbolicSrc<REAL>
{
  public:

    WxPulseTrainSrc()
            : WxHyperbolicSrc<REAL>("pulseTrain") {
    }

    void setup(const WxCryptSet& wxc) {
      // call base-class setup function
      WxHyperbolicSrc<REAL>::setup(wxc);

      // read speed of light
      _T = wxc.template get<REAL> ("Period");
      _tau = wxc.template get<REAL> ("PulseLength");
      _XC = wxc.template get<REAL> ("Xcenter");
      _YC = wxc.template get<REAL> ("Ycenter");

      _pi = 3.14159265;

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

      REAL ft = 0.0, sp;
      REAL t  = tx[0], x = tx[1], y = tx[2];
      REAL r  = sqrt(pow(x-_XC,2) + pow(y-_YC,2));
//      REAL r0 = 0.1;

      for(unsigned n = 1; n<10; ++n)
      {
          ft += 2.0/(n*_pi)*sin(_pi*n*_tau/_T)*cos(2.0*_pi*n*t/_T);
      }
      ft += _tau/_T;

      sp = ft*exp(-r*r/(0.1*0.1));
      //sp = pow(cos(_pi*r/(2.*r0)),8.);

      if(r<0.1)
          s[0] = sp;
      else
          s[0] = 0.0;

      return true;
    }

  private:
    REAL _T, _tau, _XC, _YC, _pi; // speed of light, gas_gamma
};

#endif // WXPULSETRAINSRC_H
