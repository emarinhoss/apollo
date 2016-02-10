#ifndef APCONSTANTRESISTIVITY_H
#define APCONSTANTRESISTIVITY_H

// WarpX includes
#include <wxlogger.h>

// WarpX hyperbolic solver includes
#include <wxhyperbolicsrc.h>

// std includes
#include <iostream>
#include <string>
#include <cmath>

template<class REAL>
class ApConstantResistivity : public WxHyperbolicSrc<REAL>
{
  public:

    ApConstantResistivity()
      : WxHyperbolicSrc<REAL>("constantResistivity") {
    }

    void setup(const WxCryptSet& wxc) {
      // call base-class setup function
      WxHyperbolicSrc<REAL>::setup(wxc);

      // read charge and mass of particles
      _qi = wxc.template get<REAL>("charge");
      _mi = wxc.template get<REAL>("mi");
      _me = wxc.template get<REAL>("me");

      _eta = wxc.template get<REAL>("resistivity");

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
        WxLogger::get("warpx-root.console")->
          error("*** Negative electron density in bragFrictionThermalForce src");
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
        WxLogger::get("warpx-root.console")->
          error("*** Negative electron pressure in bragFrictionThermalForce src");
        exit(1); // abort execution
      }
      if (Pe<_elcMinPres)
        Pe = _elcMinPres;
      REAL Te = Pe/(_k*ne);

      REAL rhoi  = q[5];
      if (rhoi<0.0)
      {
        WxLogger::get("warpx-root.console")->
          error("*** Negative ion density in bragFrictionThermalForce src");
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
        WxLogger::get("warpx-root.console")->
          error("*** Negative ion pressure in bragFrictionThermalForce src");
        exit(1); // abort execution
      }
      if (Pi<_ionMinPres)
        Pi = _ionMinPres;
//      REAL Ti = Pi/(_k*ni);

      REAL Rux = -_eta*ne*ne*_qe*_qe*(ue-ui);
      REAL Ruy = -_eta*ne*ne*_qe*_qe*(ve-vi);
      REAL Ruz = -_eta*ne*ne*_qe*_qe*(we-wi);

      REAL Q_delta = 3./_mi*ne*ne*_eta*_qe*_qe*(Pe/ne-Pi/ni);

      s[0] = Rux;
      s[1] = Ruy;
      s[2] = Ruz;
      s[3] = -(Rux*ui+Ruy*vi)-Q_delta;

      s[4] = -Rux;
      s[5] = -Ruy;
      s[6] = -Ruz;
      s[7] =  Q_delta;

      return true;
    }

  private:
    REAL _qi, _mi, _me, _qe, _eta;
    // ion charge, ion mass, elc mass, elc charge,
    // resistivity
    REAL _gas_gamma, _k, _elcMinPres, _ionMinPres, _pi;

};

#endif // APCONSTANTRESISTIVITY_H
