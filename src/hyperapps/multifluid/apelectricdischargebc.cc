#include "apelectricdischargebc.h"
#include <wxmath.h>
#include <algorithm>

template <typename REAL>
void
ApElectricDischargeBC<REAL>::setup(const WxCryptSet& wxc, DM dm)
{
  // call base class setup
  WxGridBC<REAL>::setup(wxc, dm);

  _Tinfty = wxc.template get<REAL>("T_infinity");
  _gammaE = wxc.template get<REAL>("secondary_emission");
  _kv = wxc.template get<REAL>("recombination_coef");
  _q = wxc.template get<REAL>("charge");
  _Efield = wxc.template get<REAL>("electric_field");
  _kB = wxc.template get<REAL>("boltzmann_constant");

  _me = wxc.template get<REAL>("electron_mass");
  _mi = wxc.template get<REAL>("ion_mass");
  _mn = wxc.template get<REAL>("neutral_mass");
  _gamma = wxc.template get<REAL>("gas_gamma");

}

template <typename REAL>
void
ApElectricDischargeBC<REAL>::applyBC(REAL *xc, REAL *nx, REAL *q, REAL *qaux, REAL *AreaInts, REAL *qBC)
{

    REAL ne = q[0]/_me;
    REAL gradTea = -2./5.*e*_Efield/_kB;

    REAL ni = q[5]/_mi;
    REAL ui = q[6]-2.*(nx[0]*q[6]+nx[1]*q[7])*nx[0]/_mi/ni;
    REAL vi = q[7]-2.*(nx[0]*q[6]+nx[1]*q[7])*nx[1]/_mi/ni;
    REAL wi = q[8]/_mi/ni;
    REAL pi = ni*_Tinfty;
    REAL Vi_perp = nx[0]*q[6]+nx[1]*q[7];
    REAL Vi_plus = std::max(Vi_perp,0.);

    REAL ex = q[10];
    REAL ey = q[11];
    REAL ez = q[12];
    REAL bx = q[13];
    REAL by = q[14];
    REAL bz = q[15];
    REAL phi= q[16];
    REAL psi= q[17];


    // electrons
    REAL Ve_par, Ve_perp = _kv-_gammaE*ni*Vi_plus/ne;

    if(Ve_perp<=0.)
        Ve_par = 0.0;
    else
        Ve_par  = (nx[0]*q[2]-nx[1]*q[1])/q[0];

    qBC[0] = q[0];
    qBC[1] = _me*ne*(Ve_perp*nx[0]-Ve_par*nx[1]);
    qBC[2] = _me*ne*(Ve_perp*nx[1]+Ve_par*nx[0]);;
    qBC[3] = q[3];
    qBC[4] = q[4];

    // ions
    qBC[5] = q[5];
    qBC[6] = ui*_mi*ni;
    qBC[7] = vi*_mi*ni;
    qBC[8] = q[8];
    qBC[9] = pi/(_gamma-1.)+0.5*ni*_mi*(ui*ui+vi*vi+wi*wi);



    // E-field
    REAL enorm = ex*nx[0] + ey*nx[1];
    REAL etang = ex*nx[1] - ey*nx[0];

    qBC[10] = enorm*nx[0] - etang*nx[1];
    qBC[11] = enorm*nx[1] + etang*nx[0];
    qBC[12] = -ez;

    // B-field
    REAL bnorm = -bx*nx[0] - by*nx[1];
    REAL btang = -bx*nx[1] + by*nx[0];

    qBC[13] = bnorm*nx[0] - btang*nx[1];
    qBC[14] = bnorm*nx[1] + btang*nx[0];
    qBC[15] = bz;

    // Divergence corrections
    qBC[16] = -phi;
    qBC[17] =  psi;
    
    // neutrals
    qBC[18] = q[18];
    qBC[19] = q[19]-2.*(nx[0]*q[19]+nx[1]*q[20])*nx[0];
    qBC[20] = q[20]-2.*(nx[0]*q[19]+nx[1]*q[20])*nx[1];
    qBC[21] = q[21];
    qBC[22] = q[22];
}

// instantiations
template class ApElectricDischargeBC<float>;
template class ApElectricDischargeBC<double>;
