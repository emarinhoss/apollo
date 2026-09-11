#ifndef APTWOFLUIDSIMPLIFIEDRMFBC_H
#define APTWOFLUIDSIMPLIFIEDRMFBC_H

// WarpX subsolver includes
#include <wxgridbc.h>

/**
 * Rotating-magnetic-field boundary condition for the two-fluid FRC problem.
 *
 * Applies, at the plasma edge r = a:
 *
 *   - a specularly reflecting wall for both fluids;
 *   - a conducting-wall condition on the in-plane electric field;
 *   - a spatially uniform transverse magnetic field of magnitude
 *     B_rmf*(1 - exp(-t/rise_time)) rotating at `frequency`, offset by `phase`;
 *   - the axial electric field that field induces, E_z = x dB_y/dt - y dB_x/dt;
 *   - an axial B_z set by flux conservation in a shell at
 *     `flux_conserver_radius`:
 *
 *         B_z(a) = [pi b^2 B_axial - Phi] / [pi (b^2 - a^2)]
 *
 *     with Phi the axial flux through the plasma, supplied by the MagneticFlux
 *     area integral. This is total-axial-flux conservation inside an ideal
 *     shell at r = b, with whatever flux the plasma no longer carries pushed
 *     into the vacuum annulus a < r < b. It assumes that annulus carries a
 *     uniform B_z, which holds because it is outside the computational domain.
 *
 * ASSUMPTIONS, and the reason this is called "simplified".
 *
 * 1. The transverse field is prescribed AT THE PLASMA EDGE. In an experiment
 *    the antenna sits outside the plasma, and the field the plasma sees is the
 *    vacuum field plus the field of the currents the plasma drives in
 *    response - the RMF is partly screened, leaving (per the TCS measurements)
 *    a mostly azimuthal field near the separatrix with a small radial
 *    component. Prescribing the vacuum field here forbids that response from
 *    reaching the edge: the antenna cannot be loaded by the plasma. The
 *    penetration threshold in B_rmf and the absorbed antenna power are
 *    therefore not the experiment's, and quantitative comparison needs the
 *    drive moved to the coils. See docs/rmf-frc-model-assessment.md.
 *
 * 2. The coil is assumed to sit outside a SLOTTED flux conserver. A solid
 *    conducting shell at these frequencies (copper skin depth ~70 um at
 *    800 kHz) screens a transverse field completely, so a coil outside a solid
 *    shell could drive nothing inside it. RMF-FRC devices segment the
 *    conserver so it passes the transverse field while still conserving axial
 *    flux; the two conditions imposed above - full B_perp, flux-conserved
 *    B_z - are exactly that arrangement.
 *
 * 3. Only the incoming characteristics of a hyperbolic system may be imposed.
 *    Supplying `c0` switches on maxwellCharacteristicGhost, which takes the
 *    outgoing ones from the interior; without it the whole prescribed field is
 *    written into the ghost state, which is over-specified. In the shipped
 *    configuration this makes no difference - the Lax-Friedrichs flux at
 *    chi = gamma = 1 is exactly upwind, and upwind discards the ghost's
 *    outgoing part - but it does away from those cleaning speeds, and in the
 *    slope limiter, which reads the ghost state directly. See
 *    src/hyperapps/maxwell/apmaxwellcharacteristics.h.
 *
 * RELATED BOUNDARY CONDITIONS, which impose DIFFERENT fields despite the
 * similar names - they are not interchangeable:
 *
 *   twoFluidRMFBC          an azimuthal field 0.5*B_rmf*(1-exp(-t/rise))
 *                          *cos(omega t + phase) * thetahat: oscillating in
 *                          time, azimuthal in space. Not a rotating uniform
 *                          transverse field, and note the factor of 0.5.
 *   twoFluidRMFHarmonicsBC the transverse field built as a harmonic series
 *                          from a coil at `coilRadius`; reads
 *                          `numberOfHarmonics` and `coilRadius`, which this
 *                          class ignores.
 *   twoFluidRMFAntennaBC   an azimuthal oscillating field, with permittivity
 *                          and permeability jumps at the wall.
 *   twoFluidOMFBC          the oscillating- rather than rotating-field
 *                          variant, with the flux-conserver term disabled.
 *
 * This class is the only one of the group that imposes the uniform transverse
 * rotating field of the classical RMF literature.
 */
template <typename REAL>
class APTwoFluidSimplifiedRMFBC : public WxGridBC<REAL>
{
  public:

/**
 * Construct a new grid-bc object
 */
    APTwoFluidSimplifiedRMFBC()
       : WxGridBC<REAL>("twoFluidSimplifiedRMFBC") {
    }

  protected:

/**
 * Setup subsolver object using supplied cryptset
 *
 * @param wxc Cryptset to use for setting
 */
    void setup(const WxCryptSet& wxc, DM dm);

/**
 * Apply BC to the upper edge along direction 'dir'
 *
 * @param dir direction in which to apply BC
 * @param arr array to which apply BC
 */
    void applyBC(REAL *xc, REAL *nx, REAL *q, REAL *qaux, REAL *AreaInts, REAL *qBC);

  private:
    REAL _omega, _baxial, _B0, _phase, _rise, _a, _b, _pi, _c0;
    bool _characteristic; // deck supplied c0; see applyBC

};

#endif // APTWOFLUIDSIMPLIFIEDRMFBC_H
