#include "aptwofluidsimplifiedrmfbc.h"
#include "maxwell/apmaxwellcharacteristics.h"
#include <wxmath.h>

template <typename REAL>
void
APTwoFluidSimplifiedRMFBC<REAL>::setup(const WxCryptSet& wxc, DM dm)
{
  // call base class setup
  WxGridBC<REAL>::setup(wxc, dm);

  REAL freq = wxc.template get<REAL>("frequency");
  _baxial = wxc.template get<REAL>("B_axial");
  _B0 = wxc.template get<REAL>("B_rmf");
  _phase = wxc.template get<REAL>("phase");
  _rise = wxc.template get<REAL>("rise_time");

  // Optional. Its presence turns on characteristic injection; see applyBC.
  _characteristic = wxc.has("c0");
  _c0 = _characteristic ? wxc.template get<REAL>("c0") : 0.0;
  _a = wxc.template get<REAL>("plasma_radius");
  _b = wxc.template get<REAL>("flux_conserver_radius");

  _pi = 3.141592653589793;

  _omega = 2*_pi*freq;


}

template <typename REAL>
void
APTwoFluidSimplifiedRMFBC<REAL>::applyBC(REAL *xc, REAL *nx, REAL *q, REAL *qaux, REAL *AreaInts, REAL *qBC)
{
    REAL x = xc[1];
    REAL y = xc[2];

    // electrons
    qBC[0] = q[0];
    qBC[1] = q[1]-2.*(nx[0]*q[1]+nx[1]*q[2])*nx[0];
    qBC[2] = q[2]-2.*(nx[0]*q[1]+nx[1]*q[2])*nx[1];
    qBC[3] = q[3];
    qBC[4] = q[4];

    // ions
    qBC[5] = q[5];
    qBC[6] = q[6]-2.*(nx[0]*q[6]+nx[1]*q[7])*nx[0];
    qBC[7] = q[7]-2.*(nx[0]*q[6]+nx[1]*q[7])*nx[1];
    qBC[8] = q[8];
    qBC[9] = q[9];

    REAL ex = q[10];
    REAL ey = q[11];

    REAL phi= q[16];
    REAL psi= q[17];

    // Area integral \int B_z \cdot dA
    // The slope limiter (tuAliabadiLimiter) calls boundary conditions to build
    // ghost states, and passes no area integrals - WxTuAliabadiLimiter::applyBc
    // hands on the NULL it was given. Every RMF boundary condition in this
    // directory read AreaInts[15] unguarded, so enabling the limiter on any
    // deck that uses one segfaulted on the first step. That is why `Limiter`
    // is commented out in the shipped rmf_frc deck.
    //
    // Falling back to the unperturbed flux pi a^2 B_axial is the right answer
    // here: the limiter only needs a ghost state to measure a slope against,
    // and the flux-conserver correction is a global quantity that does not
    // change which cells need limiting.
    REAL intBzda = AreaInts ? AreaInts[15] : 0.0;
    if(intBzda==0.0)
        intBzda = _pi*_a*_a*_baxial;

    // The field this condition wants to apply, assembled into one block so
    // that it can be handed to the characteristic injection below rather than
    // written straight into the ghost state.
    REAL qF[8];

    // E-field: conducting-wall reflection of the in-plane part.
    REAL enorm = ex*nx[0] + ey*nx[1];
    REAL etang = ex*nx[1] - ey*nx[0];

    qF[0] = enorm*nx[0] - etang*nx[1];
    qF[1] = enorm*nx[1] + etang*nx[0];

    // B-field
    //
    // The applied RMF: a spatially uniform transverse field of constant
    // magnitude B_t(t), rotating at _omega. _phase offsets the rotation.
    //
    // The phase used to be applied to the cosine component only, which does not
    // rotate the field - it changes its polarisation. At _phase = pi/2 (the
    // heavyIons deck) the two components became -B_t sin(wt) and +B_t sin(wt):
    // a linearly polarised field along a fixed axis whose magnitude swings
    // between 0 and sqrt(2) B_t, i.e. an oscillating field, not a rotating one.
    // That is a different experiment, with its own boundary condition
    // (twoFluidOMFBC). Applying the phase to both components makes it an
    // offset in the rotation, as the name says, and leaves _phase = 0 - the
    // frc2d.pin deck - bit-identical.
    REAL t = xc[0]; // current time
    REAL envelope = 1. - exp(-t/_rise);
    REAL Bt = _B0*envelope;
    REAL dBt = _B0*exp(-t/_rise)/_rise;   // d(B_t)/dt

    REAL ph = _omega*t + _phase;
    REAL Bx = -Bt*sin(ph);
    REAL By = -Bt*cos(ph);

    // Time derivatives of the applied field, for Faraday's law below.
    REAL dBx = -dBt*sin(ph) - Bt*_omega*cos(ph);
    REAL dBy = -dBt*cos(ph) + Bt*_omega*sin(ph);

    // The axial electric field induced by that rotating transverse field.
    //
    // For E = E_z zhat and a spatially uniform B_perp(t), Faraday's law
    // dB/dt = -curl E gives dE_z/dy = -dB_x/dt and dE_z/dx = +dB_y/dt, so
    //
    //     E_z = x dB_y/dt - y dB_x/dt
    //
    // which rotates with the field. The previous expression was
    // r*(...cos(wt) - B_t w sin(wt+phase)): the correct radial amplitude with
    // the azimuthal dependence dropped, so it satisfied neither component of
    // Faraday's law. E_z drives the oscillating axial currents that produce the
    // azimuthal torque, so its angular structure is the torque's structure.
    REAL Ez = x*dBy - y*dBx;

    qF[2] = Ez;
    qF[3] = Bx;
    qF[4] = By;

    REAL newBz = _b*_b*_baxial/(_b*_b-_a*_a)-intBzda/(_b*_b-_a*_a)/_pi;
//    REAL AA = -intBzda/(_b*_b-_a*_a)/_pi;
    qF[5] = newBz;
    qF[6] = -phi;
    qF[7] =  psi;

    // Impose only what a hyperbolic system allows to be imposed.
    //
    // Writing the whole of qF into the ghost state sets the outgoing
    // characteristics as well as the incoming ones, which is over-specified.
    // maxwellCharacteristicGhost keeps the outgoing ones from the interior; see
    // src/hyperapps/maxwell/apmaxwellcharacteristics.h for the decomposition
    // and test/cxx/test_maxwell_characteristics.cc for its verification.
    //
    // It makes NO DIFFERENCE at the cleaning speeds every deck in this
    // repository sets. DGnumericalFlux is Lax-Friedrichs with
    // lambda = max(c, chi c, gamma c); at chi = gamma = 1 every eigenvalue has
    // magnitude c, so that flux is exactly upwind and an upwind flux ignores
    // what the ghost says about outgoing characteristics. The two constructions
    // then agree to roundoff - measured at 2.9e-16 relative. Away from
    // chi = gamma = 1 they differ by tens of percent, and the slope limiter
    // consumes the ghost state directly rather than through a Riemann solve,
    // where the outgoing part does matter.
    //
    // Off unless the deck supplies c0, so that a deck which does not ask for
    // this is untouched down to the last bit.
    if (_characteristic)
        maxwellCharacteristicGhost(nx, _c0, q + 10, qF, qBC + 10);
    else
        for (unsigned k = 0; k < 8; ++k)
            qBC[10 + k] = qF[k];
}

// instantiations
//template class APTwoFluidSimplifiedRMFBC<float>;
template class APTwoFluidSimplifiedRMFBC<double>;
