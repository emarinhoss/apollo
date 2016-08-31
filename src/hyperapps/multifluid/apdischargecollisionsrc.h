#ifndef APDISCHARGECOLLISIONSRC_H
#define APDISCHARGECOLLISIONSRC_H

// WarpX includes
#include <wxlogger.h>

// WarpX hyperbolic solver includes
#include <wxhyperbolicsrc.h>

// std includes
#include <iostream>
#include <string>
#include <cmath>

template<class REAL>
class WxDischargeCollisionSrc : public WxHyperbolicSrc<REAL>
{
  public:

    WxDischargeCollisionSrc()
      : WxHyperbolicSrc<REAL>("electricDischargeCollisionsSrc") {
    }

    void setup(const WxCryptSet& wxc) {
      // call base-class setup function
      WxHyperbolicSrc<REAL>::setup(wxc);

      // read charge and mass of particles
      _qi = wxc.template get<REAL>("charge");
      _mi = wxc.template get<REAL>("mi");
      _me = wxc.template get<REAL>("me");
      _mn = wxc.template get<REAL>("me");
      _gamma = wxc.template get<REAL>("gas_gamma");
      _kb = wxc.template get<REAL>("boltz");

      _wi = wxc.template get<REAL>("ionization_rate");
      _wx = wxc.template get<REAL>("excitation_rate");
      _Hi = wxc.template get<REAL>("ionization_energy_loss");
      _Hx = wxc.template get<REAL>("excitation_energy_loss");

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

        // electrons
        REAL ne = q[0]/_me;
        REAL ue = q[1]/ne/_me;
        REAL ve = q[2]/ne/_me;
        REAL we = q[3]/ne/_me;
        REAL ee = q[4];
        REAL pe = (_gamma-1.)*(ee-0.5*_me*ne*(ue*ue+ve*ve+we*we));
        REAL Te = pe/ne;

        // ions
        REAL ni = q[5]/_mi;
        REAL ui = q[6]/ni/_mi;
        REAL vi = q[7]/ni/_mi;
        REAL wi = q[8]/ni/_mi;
        REAL ei = q[9];
        REAL pi = (_gamma-1.)*(ei-0.5*_mi*ni*(ui*ui+vi*vi+wi*wi));
        REAL Ti = pi/ni;

        // neutrons
        REAL nn = q[10]/_mn;
        REAL un = q[11]/nn/_mn;
        REAL vn = q[12]/nn/_mn;
        REAL wn = q[13]/nn/_mn;
        REAL en = q[14];
        REAL pn = (_gamma-1.)*(en-0.5*_mn*nn*(un*un+vn*vn+wn*wn));
        REAL Tn = pn/nn;

        // reduced masses
        REAL m_in = (_mi*_mn)/(_mi+_mn);
        REAL m_en = (_me*_mn)/(_me+_mn);

        // ion-neutral collision frequency
        REAL nu_in = 1.;

        // electron-neutral collision frequency
        REAL nu_en = 1.;

        /** Inelastic collision **/
        // Ionization
        REAL SIi = _mi*_wi;
        REAL SIe = _me*_wi;
        REAL SIn =-_mn*_wi;

        // Momentum X-fer
        // Ions
        REAL AIix = _mi*_wi*ui;
        REAL AIiy = _mi*_wi*vi;
        REAL AIiz = _mi*_wi*wi;
        // Electrons
        REAL AIex = _me*_wi*ue;
        REAL AIey = _me*_wi*ve;
        REAL AIez = _me*_wi*we;
        // Neutrals
        REAL AInx =-_mn*_wi*un;
        REAL AIny =-_mn*_wi*vn;
        REAL AInz =-_mn*_wi*wn;

        // Energy Transfer
        REAL MIi = _mi*_wi*ei;
        REAL MIe = _me*_wi*ee-_wi*_Hi-_wx*_Hx;
        REAL MIn =-_mn*_wi*en;

        /** Elastic collision **/
        // Momentum X-fer
        // Ions
        REAL AEix =-ni*m_in*nu_in*(ui-un);
        REAL AEiy =-ni*m_in*nu_in*(vi-vn);
        REAL AEiz =-ni*m_in*nu_in*(wi-wn);
        // Electrons
        REAL AEex =-ne*m_en*nu_en*(ue-un);
        REAL AEey =-ne*m_en*nu_en*(ve-vn);
        REAL AEez =-ne*m_en*nu_en*(we-wn);
        // Neutrals
        REAL AEnx =-(AEix+AEex);
        REAL AEny =-(AEiy+AEey);
        REAL AEnz =-(AEiz+AEez);

        // Energy X-fer
        REAL MEi = -(ni*m_in*nu_in)/(_mi+_mn)*(3.*(Ti-Tn)+ ( (ui-un)*(_mi*ui+_mn*un) + (vi-vn)*(_mi*vi+_mn*vn) + (wi-wn)*(_mi*wi+_mn*wn) ));
        REAL MEe = -(ne*m_en*nu_en)/(_me+_mn)*(3.*(Te-Tn)+ ( (ue-un)*(_me*ue+_mn*un) + (ve-vn)*(_me*ve+_mn*vn) + (we-wn)*(_mi*wi+_mn*wn) ));
        REAL MEn = -(MEi+MEe);

        /** Add to sources **/
        // electrons
        s[0] = SIe;
        s[1] = AIex + AEex;
        s[2] = AIey + AEey;
        s[3] = AIez + AEez;
        s[4] = MIe  + MEe;

        // ions
        s[5] = SIi;
        s[6] = AIix + AEix;
        s[7] = AIiy + AEiy;
        s[8] = AIiz + AEiz;
        s[9] = MIi  + MEi;

        // neutrals
        s[10] = SIn;
        s[11] = AInx + AEnx;
        s[12] = AIny + AEny;
        s[13] = AInz + AEnz;
        s[14] = MIn  + MEn;

      return true;
    }

  private:
    REAL _qi, _mi, _me, _mn, _qe, _gamma, _kb,
         _pi, _elcMinPres, _ionMinPres, _sgn,
         _wi, _wx, _Hi, _Hx;

};

#endif // APDISCHARGECOLLISIONSRC_H
