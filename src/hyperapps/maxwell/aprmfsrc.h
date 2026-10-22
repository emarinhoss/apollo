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

/**
 * Relaxes E_z, B_x and B_y towards a rotating transverse field outside
 * `radius`, by adding omega times the target value to their time derivatives.
 *
 * This is NOT an antenna. There is no current, so nothing here is a solution of
 * Maxwell's equations; it is a forcing term whose strength happens to be set by
 * omega, and the plasma cannot load it any more than it can load a prescribed
 * boundary field. For a drive that the plasma can screen and load, use
 * `maxwellRMFAntenna` (aprmfantennasrc.h), which sources dE_z/dt with an actual
 * current density. See docs/rmf-frc-model-assessment.md.
 *
 * Two further defects are left as they are because this class is superseded and
 * because changing them would silently change the heavyIons deck's results.
 * Both were fixed in the analogous boundary condition
 * (APTwoFluidSimplifiedRMFBC):
 *
 *   - `_phase` is applied to both components here, but as sin(wt+p) and
 *     cos(wt+p) with no minus sign, so the pair is a rotating field of the
 *     opposite handedness to the one the boundary conditions apply;
 *   - the E_z expression carries the correct radial amplitude with the
 *     azimuthal dependence dropped, and satisfies neither component of
 *     Faraday's law. The exact induced field is E_z = x dB_y/dt - y dB_x/dt.
 */
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

        // Write all three outputs unconditionally.
        //
        // WxHyperbolicSrc::compSource holds _outValues as a member, does not
        // clear it between calls, and adds it to the state:
        //
        //     sfull[_outIndices[i]] += _outValues[i];
        //
        // so the `if (r > _r0)` below used to leave the PREVIOUS quadrature
        // point's values in place, and they were added at every point inside
        // r0 - the whole plasma interior of the heavyIons deck, which is the
        // only deck that uses this class (r0 = 0.025 against a plasma radius
        // of 0.030). Whatever the last point outside r0 happened to produce
        // was injected into E_z, B_x and B_y across the column.
        s[0] = 0.0;
        s[1] = 0.0;
        s[2] = 0.0;

        if (r <= _r0)
            return true;

        // RMF
        REAL t = tx[0]; // current time
        REAL Bt = _B0*(1.-exp(-t/_rise));

        REAL Br = Bt*sin(_omega*t+_phase);
        REAL Bc = Bt*cos(_omega*t+_phase);

        REAL Ez = r*(_B0*(-exp(-t/_rise))/_rise*cos(_omega*t+_phase)-Bt*_omega*sin(_omega*t+_phase));

        s[1] = _omega*Br;
        s[2] = _omega*Bc;
        s[0] = _omega*Ez;

      return true;
    }

  private:
    REAL _omega, _B0, _pi, _rise, _phase, _r0;

};

#endif // APRMFSRC_H
