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

      // set minimum electron pressure to prevent negative
      if (wxc.has("elcMinDensity"))
        _elcMinDens = wxc.template get<REAL>("elcMinDensity");
      else
        _elcMinDens = 0.0;

      // set minimum ion pressure to prevent negative
      if (wxc.has("ionMinDensity"))
        _ionMinDens = wxc.template get<REAL>("ionMinDensity");
      else
        _ionMinDens = 0.0;


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
      if (rhoe<_elcMinDens && _elcMinDens==0.0)
      {
        WxLogger::get("apollo-root.console")->
          error("*** Negative electron density in ConstantResistivity src ***\n");
        exit(1); // abort execution
      }
      if(rhoe<_elcMinDens)
          rhoe = _elcMinDens;

      REAL ne = rhoe/_me;
      REAL ue = q[1]/rhoe;
      REAL ve = q[2]/rhoe;
      REAL we = q[3]/rhoe;
      REAL Ee = q[4];
      REAL Pe = (_gas_gamma-1)*(Ee-0.5*rhoe*(ue*ue+ve*ve+we*we));
      if (Pe<_elcMinPres && _elcMinPres==0.0)
      {
        WxLogger::get("apollo-root.console")->
          error("*** Negative electron pressure in ConstantResistivity src ***\n");
        exit(1); // abort execution
      }
      if (Pe<_elcMinPres)
        Pe = _elcMinPres;
//      REAL Te = Pe/(_k*ne);

      REAL rhoi  = q[5];
      if (rhoi<_ionMinDens && _ionMinDens==0.0)
      {
        WxLogger::get("apollo-root.console")->
          error("*** Negative ion density in ConstantResistivity src ***\n");
        exit(1); // abort execution
      }
      if(rhoi<_ionMinDens)
          rhoi = _ionMinDens;

      REAL ni = rhoi/_mi;
      REAL ui = q[6]/rhoi;
      REAL vi = q[7]/rhoi;
      REAL wi = q[8]/rhoi;
      REAL Ei = q[9];
      REAL Pi = (_gas_gamma-1)*(Ei-0.5*rhoi*(ui*ui+vi*vi+wi*wi));
      if (Pi<_ionMinPres && _ionMinPres==0.0)
      {
        WxLogger::get("apollo-root.console")->
          error("*** Negative ion pressure in ConstantResistivity src ***\n");
        exit(1); // abort execution
      }
      if (Pi<_ionMinPres)
        Pi = _ionMinPres;
//      REAL Ti = Pi/(_k*ni);

      REAL Rux = -_eta*ne*ne*_qe*_qe*(ue-ui);
      REAL Ruy = -_eta*ne*ne*_qe*_qe*(ve-vi);
      REAL Ruz = -_eta*ne*ne*_qe*_qe*(we-wi);

      REAL Q_delta = 3./_mi*ne*ne*_eta*_qe*_qe*(Pe/ne-Pi/ni);

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
      // Friction is three-dimensional here, so the work term must be too:
      // omitting Ruz*wi would break conservation even with the signs right.
      REAL Rdotui = Rux*ui + Ruy*vi + Ruz*wi;

      s[0] = Rux;
      s[1] = Ruy;
      s[2] = Ruz;
      s[3] =  Rdotui - Q_delta;

      s[4] = -Rux;
      s[5] = -Ruy;
      s[6] = -Ruz;
      s[7] = -Rdotui + Q_delta;

      return true;
    }

  private:
    REAL _qi, _mi, _me, _qe, _eta;
    // ion charge, ion mass, elc mass, elc charge,
    // resistivity
    REAL _gas_gamma, _k, _elcMinPres, _ionMinPres, _elcMinDens, _ionMinDens, _pi;

};

#endif // APCONSTANTRESISTIVITY_H
