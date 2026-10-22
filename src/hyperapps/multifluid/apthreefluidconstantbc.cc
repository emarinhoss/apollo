#include "apthreefluidconstantbc.h"

template <typename REAL>
void
ApThreeFluidConstantBC<REAL>::setup(const WxCryptSet& wxc, DM dm)
{
    // call base class setup
    WxGridBC<REAL>::setup(wxc, dm);

    // set left initial state
    _pe = wxc.template get<REAL>("electron_pressure");
    _pi = wxc.template get<REAL>("ion_pressure");
    _pn = wxc.template get<REAL>("neutral_pressure");

    _ne = wxc.template get<REAL>("electron_numDens");
    _ni = wxc.template get<REAL>("ion_numDens");
    _nn = wxc.template get<REAL>("neutral_numDens");

    _me = wxc.template get<REAL>("electron_mass");
    _mi = wxc.template get<REAL>("ion_mass");
    _mn = wxc.template get<REAL>("neutral_mass");

    _gas_gamma = wxc.template get<REAL>("gas_gamma");

}

template <typename REAL>
void
ApThreeFluidConstantBC<REAL>::applyBC(REAL *xc, REAL *nx, REAL *q, REAL *qaux, REAL *AreaInts, REAL *qBC)
{
    // electrons
    qBC[0] = _me*_ne;
    qBC[1] = 0.;
    qBC[2] = 0.;
    qBC[3] = 0.;
    qBC[4] = _pe/(_gas_gamma-1.);

    // ions
    qBC[5] = _mi*_ni;
    qBC[6] = 0.;
    qBC[7] = 0.;
    qBC[8] = 0.;
    qBC[9] = _pi/(_gas_gamma-1.);

    // neutrals
    qBC[10] = _mn*_nn;
    qBC[11] = 0.;
    qBC[12] = 0.;
    qBC[13] = 0.;
    qBC[14] = _pn/(_gas_gamma-1.);


}

// instantiations
//template class ApThreeFluidConstantBC<float>;
template class ApThreeFluidConstantBC<double>;
