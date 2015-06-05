#ifndef __wxlorentzforcesrc__h__
#define __wxlorentzforcesrc__h__

// WarpX includes

// WarpX hyperbolic solver includes
#include <wxhyperbolicsrc.h>

// std includes
#include <iostream>
#include <string>

template<class REAL>
class WxLorentzForceSrc : public WxHyperbolicSrc<REAL>
{
  public:

    WxLorentzForceSrc()
        : WxHyperbolicSrc<REAL>("lorentzForces") {
    }

    void setup(const WxCryptSet& wxc) 
    {
      // call base-class setup function
      WxHyperbolicSrc<REAL>::setup(wxc);
      // read charge and mass of particles
      _q = wxc.template get<REAL>("charge");
      _m = wxc.template get<REAL> ("mass");
      _qbym = _q/_m;
    }
    
  /**
   * Calculate value of Lorentz force source terms at given time/space
   * locations for given input conserved variables.
   *
   * @param n Length of q array
   * @param tx[0] is time and tx[1..3] are spatial location
   * @param q Conserved variables at tx
   * @param qaux Conserved variables at tx
   * @param s Source term
   */
    bool src(unsigned n, REAL *tx, REAL *q, REAL *qaux, REAL *s) 
    {
      // assumes that q is [rho, rho*u, rho*v, rho*w, Ex, Ey, Ez, Bx, By, Bz]

      REAL rho  = q[0];
      REAL rhou = q[1];
      REAL rhov = q[2];
      REAL rhow = q[3];
      REAL ex = q[4];
      REAL ey = q[5];
      REAL ez = q[6];
      REAL bx = q[7];
      REAL by = q[8];
      REAL bz = q[9];

      // compute source terms
      s[0] = _qbym*(rho*ex+rhov*bz-rhow*by); // x-momentum
      s[1] = _qbym*(rho*ey+rhow*bx-rhou*bz); // y-momentum
      s[2] = _qbym*(rho*ez+rhou*by-rhov*bx); // z-momentum
      s[3] = _qbym*(rhou*ex+rhov*ey+rhow*ez); // energy

      return true;
    }

  private:
    REAL _q, _m, _qbym; // charge, mass and charge/mass ratio
};

#endif //  __wxlorentzforcesrc__h__
