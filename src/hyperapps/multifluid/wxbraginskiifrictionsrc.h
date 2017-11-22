#ifndef WXBRAGINSKIIFRICTIONSRC_H
#define WXBRAGINSKIIFRICTIONSRC_H

// WarpX includes
#include <wxlogger.h>

// WarpX hyperbolic solver includes
#include <wxhyperbolicsrc.h>

// std includes
#include <iostream>
#include <string>
#include <cmath>

template<class REAL>
class WxBragFrictionSrc : public WxHyperbolicSrc<REAL>
{
  public:

    WxBragFrictionSrc()
      : WxHyperbolicSrc<REAL>("bragFriction") {
    }

    void setup(const WxCryptSet& wxc) {
      // call base-class setup function
      WxHyperbolicSrc<REAL>::setup(wxc);

      // read charge and mass of particles
      _qi = wxc.template get<REAL>("charge");
      _mi = wxc.template get<REAL>("mi");
      _me = wxc.template get<REAL>("me");
      _eps0 = wxc.template get<REAL>("epsilon0");
      _gas_gamma = wxc.template get<REAL>("gas_gamma");
      _k = wxc.template get<REAL>("boltz");

      // set minimum electron pressure to prevent negative
      if (wxc.has("elcMinPressure"))
        _elcMinPres = wxc.template get<REAL>("elcMinPressure");
      else
        _elcMinPres = 0.0;

      // set minimum ion pressure to prevent negative
      if (wxc.has("ionMinPressure"))
        _ionMinPres = wxc.template get<REAL>("ionMinPressure");
      else
        _ionMinPres = 0.0;


      _qe = -1.0*_qi;
      _qbymi = _qi/_mi;
      _qbyme = _qe/_me;

      _pi = 3.14159265;
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
      // assumes that q is conserved two-fluid variables:
      // [rhoe, rhoe*ue, rhoe*ve, rhoe*we, Ee,
      //  rhoi, rhoi*ui, rhoi*vi, rhoi*wi, Ei,
      //  Bx, By, Bz]

      // Source terms here: momentum transport
      // Rei from Braginskii

      // *** ALL VALUES FOR CONSTANTS ASSUME Z=1 ***

      REAL rhoe  = q[0];
      if (rhoe<0.0)
      {
        WxLogger::get("apollo-root.console")->
          error("*** Negative electron density in bragFriction src");
        exit(1); // abort execution
      }
      REAL ne = rhoe/_me;
      REAL ue = q[1]/rhoe;
      REAL ve = q[2]/rhoe;
      REAL we = q[3]/rhoe;
      REAL Ee = q[4];
      REAL Pe = (_gas_gamma-1)*(Ee-0.5*rhoe*(ue*ue+ve*ve+we*we));
      if (Pe<_elcMinPres && _elcMinPres==0.0)
      {
        WxLogger::get("apollo-root.console")->
          error("*** Negative electron pressure in bragFriction src");
        exit(1); // abort execution
      }
      if (Pe<_elcMinPres)
        Pe = _elcMinPres;
      REAL Te = Pe/(_k*ne);
      REAL Vte= sqrt(2.*_k*Te/_me);

      REAL rhoi  = q[5];
      if (rhoi<0.0)
      {
        WxLogger::get("apollo-root.console")->
          error("*** Negative ion density in bragFriction src");
        exit(1); // abort execution
      }
      REAL ni = rhoi/_mi;
      REAL ui = q[6]/rhoi;
      REAL vi = q[7]/rhoi;
      REAL wi = q[8]/rhoi;
      REAL Ei = q[9];
      REAL Pi = (_gas_gamma-1)*(Ei-0.5*rhoi*(ui*ui+vi*vi+wi*wi));
      if (Pi<_ionMinPres && _ionMinPres==0.0)
      {
        WxLogger::get("apollo-root.console")->
          error("*** Negative ion pressure in bragFriction src");
        exit(1); // abort execution
      }
      if (Pi<_ionMinPres)
        Pi = _ionMinPres;

      REAL u = ue-ui;
      REAL v = ve-vi;
      REAL w = we-wi;
//      REAL vi_mag = sqrt(ui*ui+vi*vi+wi*wi);

      REAL bx = q[10];
      REAL by = q[11];
      REAL bz = q[12];
      REAL bmag = sqrt(bx*bx + by*by + bz*bz);

      // compute electron-electron collision frequency
      REAL edebye = sqrt(_eps0*_k*Te/(ne*_qe*_qe));
      REAL loglambda = log(4./3.*_pi*ne*pow(edebye,3.0));
//      REAL nue = (pow(_qe,4.0)*ne*loglambda)/(4.*_pi*_eps0*_eps0*_mi*_me)/(pow(vi_mag,3)+1.3*pow(Vte,3));
      REAL nue = 8.06e5*ni*loglambda/pow(Vte,3);

      if (nue<0)
      {
        WxLogger::get("apollo-root.console")->
          error("*** Negative electron collision frequency in bragFriction src");
        exit(1); // abort execution
      }

      // compute electron cyclotron frequency
      REAL wce = -_qe*bmag/_me;

      // Compute variables used in following formulations
      // same variable names used as in Braginskii paper
      REAL xe = wce/nue;

      REAL delta0 = 3.7703;
      REAL delta1 = 14.79;
      REAL deltae = pow(xe,4.0) + delta1*xe*xe + delta0;

      // ----------------------------------------------
      // ------- Momentum transport source term -------
      // ----------------------------------------------
      // Rei = Ru -> frictional force

      // Step 1: frictional force term, Ru
      // compute coefficients
      REAL Rux, Ruy, Ruz;
      if (bmag==0)
      {
        REAL alpha0 = 0.5129;
        REAL alpha_par = rhoe*nue*alpha0;

        Rux = -alpha_par*u;
        Ruy = -alpha_par*v;
        Ruz = -alpha_par*w;

      }
      else
      {
        REAL alpha0 = 0.5129;
        REAL alpha0p = 1.837;
        REAL alpha1p = 6.416;
        REAL alpha0pp = 0.7796;
        REAL alpha1pp = 1.704;

        REAL alpha_par = rhoe*nue*alpha0;
        REAL alpha_perp = rhoe*nue*(1.0-(alpha1p*xe*xe+alpha0p)/deltae);
        REAL alpha_cross = rhoe*nue*xe*(alpha1pp*xe*xe+alpha0pp)/deltae;

        // compute relative velocity parallel to magnetic field
        REAL upar = (u*bx/bmag + v*by/bmag + w*bz/bmag)*bx/bmag;
        REAL vpar = (u*bx/bmag + v*by/bmag + w*bz/bmag)*by/bmag;
        REAL wpar = (u*bx/bmag + v*by/bmag + w*bz/bmag)*bz/bmag;

        // compute relative velocity perpendicular to magnetic field
        REAL uperp = u-upar;
        REAL vperp = v-vpar;
        REAL wperp = w-wpar;

        // compute relative velocity in direction perpendicular to
        // both magnetic field and velocity
        REAL ubcrossu = by/bmag*w - bz/bmag*v;
        REAL vbcrossu = bz/bmag*u - bx/bmag*w;
        REAL wbcrossu = bx/bmag*v - by/bmag*u;

        Rux = -alpha_par*upar - alpha_perp*uperp + alpha_cross*ubcrossu;
        Ruy = -alpha_par*vpar - alpha_perp*vperp + alpha_cross*vbcrossu;
        Ruz = -alpha_par*wpar - alpha_perp*wperp + alpha_cross*wbcrossu;
      }

      REAL Q_delta = 3*_mi/_me*ne*nue*(Pe/ne-Pi/ni);

      s[0] = Rux;
      s[1] = Ruy;
      s[2] = Ruz;
      s[3] = -(Rux*ui+Ruy*vi+Ruz*wi)-Q_delta;

      s[4] = -Rux;
      s[5] = -Ruy;
      s[6] = -Ruz;
      s[7] = Q_delta;

      return true;
    }

  private:
    REAL _qi, _mi, _me, _qe, _qbymi, _qbyme,
        _eps0, _gas_gamma, _k, _pi,
        _elcMinPres, _ionMinPres, _sgn;
    bool _norm;
    // ion charge, ion mass, elc charge, elc mass,
    // ion charge-to-mass, elc charge-to-mass,
    // epsilon0, boltzmann const, gas gamma, pi,
    // elc min pressure, ion min pressure

};

#endif // WXBRAGINSKIIFRICTIONSRC_H
