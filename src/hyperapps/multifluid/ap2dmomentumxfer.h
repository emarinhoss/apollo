#ifndef AP2DMOMENTUMXFER_H
#define AP2DMOMENTUMXFER_H

// WarpX includes
#include <wxlogger.h>

// WarpX hyperbolic solver includes
#include <wxhyperbolicsrc.h>

// std includes
#include <iostream>
#include <string>
#include <cmath>

template<class REAL>
class Ap2DMomentumXfer : public WxHyperbolicSrc<REAL>
{
  public:

    Ap2DMomentumXfer()
      : WxHyperbolicSrc<REAL>("MomentumXfer2D") {
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
          error("*** Negative electron density in 2dMomentumXfer src");
        exit(1); // abort execution
      }
      REAL ne = rhoe/_me;
      REAL ue = q[1]/rhoe;
      REAL ve = q[2]/rhoe;
      REAL Ee = q[4];
      REAL Pe = (_gas_gamma-1)*(Ee-0.5*rhoe*(ue*ue+ve*ve));
      if (Pe<_elcMinPres && _elcMinPres==0.0)
      {
        WxLogger::get("apollo-root.console")->
          error("*** Negative electron pressure in 2dMomentumXfer src");
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
          error("*** Negative ion density in 2dMomentumXfer src");
        exit(1); // abort execution
      }
      REAL ni = rhoi/_mi;
      REAL ui = q[6]/rhoi;
      REAL vi = q[7]/rhoi;
      REAL Ei = q[9];
      REAL Pi = (_gas_gamma-1)*(Ei-0.5*rhoi*(ui*ui+vi*vi));
      if (Pi<_ionMinPres && _ionMinPres==0.0)
      {
        WxLogger::get("apollo-root.console")->
          error("*** Negative ion pressure in 2dMomentumXfer src");
        exit(1); // abort execution
      }
      if (Pi<_ionMinPres)
        Pi = _ionMinPres;

      REAL u = ue-ui;
      REAL v = ve-vi;

      REAL bx = q[10];
      REAL by = q[11];
      REAL bmag = sqrt(bx*bx + by*by);

      // compute electron-electron collision frequency
      REAL edebye = sqrt(_eps0*_k*Te/(ne*_qe*_qe));
      REAL loglambda = log(4./3.*_pi*ne*pow(edebye,3.0));
//      REAL nue = (pow(_qe,4.0)*ne*loglambda)/(4.*_pi*_eps0*_eps0*_mi*_me)/(pow(vi_mag,3)+1.3*pow(Vte,3));
      REAL nue = 8.06e5*ni*loglambda/pow(Vte,3);

      if (nue<0)
      {
        WxLogger::get("apollo-root.console")->
          error("*** Negative electron collision frequency in 2dMomentumXfer src");
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
      REAL Rux, Ruy;
      if (bmag==0)
      {
        REAL alpha0 = 0.5129;
        REAL alpha_par = rhoe*nue*alpha0;

        Rux = -alpha_par*u;
        Ruy = -alpha_par*v;

      }
      else
      {
        REAL alpha0 = 0.5129;
        REAL alpha0p = 1.837;
        REAL alpha1p = 6.416;

        REAL alpha_par = rhoe*nue*alpha0;
        REAL alpha_perp = rhoe*nue*(1.0-(alpha1p*xe*xe+alpha0p)/deltae);

        // compute relative velocity parallel to magnetic field
        REAL upar = (u*bx/bmag + v*by/bmag)*bx/bmag;
        REAL vpar = (u*bx/bmag + v*by/bmag)*by/bmag;

        // compute relative velocity perpendicular to magnetic field
        REAL uperp = u-upar;
        REAL vperp = v-vpar;

        // compute relative velocity in direction perpendicular to
        // both magnetic field and velocity
        Rux = -alpha_par*upar - alpha_perp*uperp;
        Ruy = -alpha_par*vpar - alpha_perp*vperp;
      }

      // Braginskii electron-ion thermal equilibration,
      //     Q_Delta = 3 (m_e/m_i) n_e nu_ei (T_e - T_i).
      // The mass ratio was inverted here (3*_mi/_me), overstating the rate by
      // (m_i/m_e)^2 - a factor of 3.4e6 for hydrogen, which would equilibrate
      // the species essentially instantaneously instead of making it the
      // slowest collisional process in the system. The same coefficient in
      // apconstantresistivity.h, written as (3/m_i) n^2 eta e^2 with
      // eta = m_e nu_ei/(n e^2), reduces to this and has always been right.
      REAL Q_delta = 3*_me/_mi*ne*nue*(Pe/ne-Pi/ni);

      // Collisional energy exchange, in total-energy variables.
      //
      // With R the friction on the electrons (so -R on the ions) and
      // w = u_e - u_i, the sources are s_e = R.u_e + Q_e and
      // s_i = -R.u_i + Q_i, and conservation forces Q_e + Q_i = -R.w, the
      // frictional heat (equal to eta*J^2 for the isotropic law). Braginskii
      // deposits that in the electrons and moves Q_Delta from electrons to
      // ions, so Q_e = -R.w - Q_Delta and Q_i = Q_Delta. Substituting, both
      // collapse to one dot product with the ION velocity:
      //
      //     s_e = +R.u_i - Q_Delta ,   s_i = -R.u_i + Q_Delta
      //
      // which sum to zero for any friction law. The frictional heating needs
      // no separate term: the momentum source already removes exactly that
      // much kinetic energy from the electron fluid, and leaving it out of the
      // energy source is what turns it into heat.
      //
      // Previously s_e carried -R.u_i and s_i carried no work term at all,
      // giving total energy a spurious source of -R.u_i.
      // Only the in-plane friction is applied here, so only in-plane work.
      REAL Rdotui = Rux*ui + Ruy*vi;

      s[0] = Rux;
      s[1] = Ruy;
      s[2] = 0.0;
      s[3] =  Rdotui - Q_delta;

      s[4] = -Rux;
      s[5] = -Ruy;
      s[6] =  0.0;
      s[7] = -Rdotui + Q_delta;

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

#endif // AP2DMOMENTUMXFER_H
