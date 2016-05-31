#ifndef APRMFELECTRONVELOCITYSRC_H
#define APRMFELECTRONVELOCITYSRC_H

// WarpX includes
#include <wxlogger.h>

// WarpX hyperbolic solver includes
#include <wxhyperbolicsrc.h>

// std includes
#include <iostream>
#include <string>
#include <cmath>

template<class REAL>
class ApRMFelectronVelocitySrc : public WxHyperbolicSrc<REAL>
{
  public:

    ApRMFelectronVelocitySrc()
      : WxHyperbolicSrc<REAL>("rmfElectronVelSrc") {
    }

    void setup(const WxCryptSet& wxc) {
      // call base-class setup function
      WxHyperbolicSrc<REAL>::setup(wxc);

      // read charge and mass of particles
      REAL freq = wxc.template get<REAL>("frequency");
      _r0 = wxc.template get<REAL>("radius");

      REAL pi = 3.141592653589793;

      _omega = 2.*pi*freq;
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

        //
        REAL rho = q[0];
//        REAL Vx  = q[1]/rho;
//        REAL Vy  = q[2]/rho;

//        REAL Vr = Vx*cos(theta) + Vy*sin(theta);
//        REAL Vt =-Vx*sin(theta) + Vy*cos(theta);

        REAL Vel = _omega*r;

        if(r>_r0){
            s[0] =-Vel*sin(theta)*rho;
            s[1] = Vel*cos(theta)*rho;
        }

      return true;
    }

  private:
    REAL _omega, _r0;

};

#endif // APRMFELECTRONVELOCITYSRC_H
