#ifndef APMAGNETICFLUXCALC_H
#define APMAGNETICFLUXCALC_H

// Apollo hyperbolic solver includes
#include <apareaintegral.h>

// std includes
#include <iostream>
#include <string>

template<class REAL>
class ApMagneticFluxCalc : public ApAreaIntegral<REAL>
{
public:

    ApMagneticFluxCalc()
            : ApAreaIntegral<REAL>("MagneticFlux")
      {
      }

    void setup(const WxCryptSet& wxc) {
        // call base-class setup function
        ApAreaIntegral<REAL>::setup(wxc);

        _pi = 3.141592653589793;
        REAL radius = wxc.template get<REAL>("plasma_radius");

        _area = _pi*radius*radius;

    }

/**
 * Calculate value of a surface area term at given time/space locations for
 * given input conserved variables.
 *
 * @param n Length of q array
 * @param tx[0] is time and tx[1..3] are spatial location
 * @param q Conserved variables at tx
 * @param qaux Auxillary variables
 * @param s area integral term
 */
    bool integral(unsigned n, REAL *tx, REAL *q, REAL *qaux, REAL *s) {
        // assumes that q is [rho*u, rho*v, rho*w]

        s[0] = q[0];

        return true;
    }

private:
    REAL _pi, _area;

};

#endif // APMAGNETICFLUXCALC_H
