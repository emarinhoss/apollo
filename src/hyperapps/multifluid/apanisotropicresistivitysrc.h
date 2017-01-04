#ifndef APANISOTROPICRESISTIVITYSRC_H
#define APANISOTROPICRESISTIVITYSRC_H

// WarpX includes
#include <wxlogger.h>

// WarpX hyperbolic solver includes
#include <wxhyperbolicsrc.h>

// std includes
#include <iostream>
#include <string>
#include <cmath>

template<class REAL>
class ApAnisotropicResistivitySrc : public WxHyperbolicSrc<REAL>
{
  public:

    ApAnisotropicResistivitySrc()
      : WxHyperbolicSrc<REAL>("anisotropicResistivity") {
    }

    void setup(const WxCryptSet& wxc) {
      // call base-class setup function
      WxHyperbolicSrc<REAL>::setup(wxc);

      // read charge and mass of particles
      _qi = wxc.template get<REAL>("charge");
      _mi = wxc.template get<REAL>("mi");
      _me = wxc.template get<REAL>("me");

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


      _pi = 3.14159265;

      _qe = -1.0*_qi;
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
      //  rhoi, rhoi*ui, rhoi*vi, rhoi*wi, Ei]

      REAL rhoe  = q[0];
      if (rhoe<0.0)
      {
        WxLogger::get("apollo-root.console")->
          error("*** Negative electron density in AnisotropicResistivity src");
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
          error("*** Negative electron pressure in AnisotropicResistivity src");
        exit(1); // abort execution
      }
      if (Pe<_elcMinPres)
        Pe = _elcMinPres;
      REAL Te = Pe/ne/_k/11600.;

      REAL rhoi  = q[5];
      if (rhoi<0.0)
      {
        WxLogger::get("apollo-root.console")->
          error("*** Negative ion density in AnisotropicResistivity src");
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
          error("*** Negative ion pressure in AnisotropicResistivity src");
        exit(1); // abort execution
      }
      if (Pi<_ionMinPres)
        Pi = _ionMinPres;

      // Magnetic field
      REAL bx = q[10];
      REAL by = q[11];
//      REAL bz = q[12];
      REAL bmag = sqrt(bx*bx + by*by);

      REAL loglambda = 10.;
      REAL eta_par = 5.2e-5*loglambda/pow(Te,1.5);
      REAL eta_perp= 2.0*eta_par;

      // velocity differences
      REAL u = ue-ui;
      REAL v = ve-vi;

      REAL Rux, Ruy;

      if(bmag == 0.0)
      {
          Rux = -eta_par*ne*ne*_qe*_qe*u;
          Ruy = -eta_par*ne*ne*_qe*_qe*v;
      }
      else
      {
          // compute relative velocity parallel to magnetic field
          REAL upar = (u*bx/bmag + v*by/bmag)*bx/bmag;
          REAL vpar = (u*bx/bmag + v*by/bmag)*by/bmag;

          // compute relative velocity perpendicular to magnetic field
          REAL uperp = u-upar;
          REAL vperp = v-vpar;

          REAL Rpar_x =  -eta_par*ne*ne*_qe*_qe*upar;
          REAL Rpar_y =  -eta_par*ne*ne*_qe*_qe*vpar;

          REAL Rperp_x = -eta_perp*ne*ne*_qe*_qe*uperp;
          REAL Rperp_y = -eta_perp*ne*ne*_qe*_qe*vperp;

          Rux = Rpar_x + Rperp_x;
          Ruy = Rpar_y + Rperp_y;
      }

      REAL Q_delta = 3./_mi*ne*ne*eta_par*_qe*_qe*(Pe/ne-Pi/ni);

      s[0] = Rux;
      s[1] = Ruy;
      s[2] = 0.0;
      s[3] = -(Rux*ui+Ruy*vi)-Q_delta;

      s[4] = -Rux;
      s[5] = -Ruy;
      s[6] =  0.0;
      s[7] =  Q_delta;

      return true;
    }

  private:
    REAL _qi, _mi, _me, _qe;
    // ion charge, ion mass, elc mass, elc charge,
    // resistivity
    REAL _gas_gamma, _k, _elcMinPres, _ionMinPres, _pi;

};

#endif // APANISOTROPICRESISTIVITYSRC_H
