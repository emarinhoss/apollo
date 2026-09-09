#ifndef APMAXWELLCHARACTERISTICS_H
#define APMAXWELLCHARACTERISTICS_H

/**
 * Characteristic-consistent boundary states for the perfectly-hyperbolic
 * Maxwell system. Phase 2 of docs/rmf-frc-model-assessment.md.
 *
 * THE PROBLEM. A boundary condition that prescribes a field writes the whole
 * of it into the ghost state - all six field components plus the two cleaning
 * potentials. That over-specifies a hyperbolic system. Only the characteristics
 * travelling INTO the domain may be imposed; the ones travelling out are
 * determined by the interior and imposing them reflects outgoing waves back in.
 * Every RMF boundary condition in this tree did exactly that.
 *
 * THE DECOMPOSITION. With the conserved state ordered
 * [E_x, E_y, E_z, B_x, B_y, B_z, phi, psi] and the flux of WxPHMaxwellEqn,
 *
 *     f_x = [chi c^2 phi, c^2 B_z, -c^2 B_y, gamma psi, -E_z, E_y,
 *            chi E_x, gamma c^2 B_x]
 *
 * the flux Jacobian in a frame whose x-axis is the outward face normal n has
 * eigenvalues +-c (twice each), +-chi c and +-gamma c, with left eigenvectors
 * that are pairwise sums and differences. Writing the field in the face frame
 * (n; t1 = n rotated a quarter turn in the plane; t2 = z):
 *
 *     OUTGOING (speed > 0, leaving)   INCOMING (speed < 0, entering)
 *       v1 = E_t1 + c B_t2   (+c)       w1 = E_t1 - c B_t2   (-c)
 *       v2 = E_t2 - c B_t1   (+c)       w2 = E_t2 + c B_t1   (-c)
 *       v3 = E_n   + c phi   (+chi c)   w3 = E_n   - c phi   (-chi c)
 *       v4 = B_n   + psi/c   (+gam c)   w4 = B_n   - psi/c   (-gam c)
 *
 * All eight were checked against the Jacobian assembled numerically from
 * WxPHMaxwellEqn::flux itself; test/cxx/test_maxwell_characteristics.cc repeats
 * that check, so a change to the flux cannot silently invalidate this.
 *
 * Note that chi and gamma do not appear below. They scale the two cleaning wave
 * speeds but not their signs, and it is only the sign that decides which
 * characteristics are incoming - so the ghost state is the same for any
 * positive pair.
 *
 * THE GHOST STATE. Take the outgoing invariants from the interior and the
 * incoming ones from the prescribed field. Equivalently, and this is how it is
 * computed below because it is far better conditioned, add to the interior
 * state the incoming part of (prescribed - interior). An exact Riemann solve
 * at the face then sees precisely the intended boundary state, and:
 *
 *   - when the prescribed field equals the interior state the ghost IS the
 *     interior state, i.e. this degenerates to the zero-gradient outflow
 *     condition (phmOpenBC), as a boundary injection should;
 *   - when the prescribed field carries only incoming characteristics it is
 *     passed through exactly.
 *
 * A NOTE ON THE FLUX. WxPHMaxwellEqn::DGnumericalFlux is Lax-Friedrichs with
 * lambda = max(c, chi c, gamma c). When chi = gamma = 1 - what every deck in
 * this repository sets - every eigenvalue has magnitude c, so |A| = c I and
 * that Lax-Friedrichs flux IS the exact upwind flux. The construction below is
 * then not merely consistent but exact. For chi or gamma other than 1 the flux
 * is more diffusive than upwind and the realised boundary is correspondingly
 * smeared; the ghost state is still the right one to supply.
 *
 * WHY NOT eigenSystem(). The plan proposed injecting through
 * WxPHMaxwellEqn::eigenSystem, which assembles left and right eigenvectors.
 * That function is never called: the only caller is
 * WxHyperbolicEqnSet::eigenSystem, which nothing calls in turn. It is also
 * wrong where it is exercised by eye - it reads the cleaning speeds as
 * `gamma = q[6]; kappa = q[7]`, which are the cleaning POTENTIALS phi and psi,
 * not the speeds - and it is axis-aligned, so an unstructured face normal would
 * need a rotation around it anyway. The closed form above is shorter than the
 * rotation would be, and is tested.
 */

/**
 * Ghost state for a prescribed-field boundary, imposing only what may be
 * imposed.
 *
 * @param nx      outward unit normal of the face, (n_x, n_y)
 * @param c0      speed of light as the equation object uses it
 * @param qInt    interior Maxwell state, 8 components from E_x
 * @param qBnd    the field the boundary condition wants to apply, same layout
 * @param qGhost  out: the ghost state to hand to the numerical flux. May alias
 *                neither qInt nor qBnd.
 */
template <typename REAL>
inline void
maxwellCharacteristicGhost(const REAL *nx, REAL c0,
                           const REAL *qInt, const REAL *qBnd, REAL *qGhost)
{
    const REAL n0 = nx[0], n1 = nx[1];

    // Work in the DIFFERENCE between the prescribed field and the interior,
    // not in the invariants themselves. Both forms are algebraically the same;
    // this one is well conditioned and the other is not. Reconstructing psi
    // from (v4 - w4) c/2 subtracts two numbers that are each about B_n and
    // differ by 2 psi/c, so with the reduced speed of light and psi small it
    // loses seven digits - enough that feeding the interior state back in did
    // not return it. Here a zero difference gives exactly zero correction.
    const REAL dEx = qBnd[0] - qInt[0], dEy = qBnd[1] - qInt[1];
    const REAL dBx = qBnd[3] - qInt[3], dBy = qBnd[4] - qInt[4];

    // Into the face frame. t1 = (-n1, n0) is n turned a quarter turn
    // anticlockwise in the plane; t2 is z, which is already a component.
    const REAL dEn  =  dEx*n0 + dEy*n1;
    const REAL dEt1 = -dEx*n1 + dEy*n0;
    const REAL dEt2 =  qBnd[2] - qInt[2];
    const REAL dBn  =  dBx*n0 + dBy*n1;
    const REAL dBt1 = -dBx*n1 + dBy*n0;
    const REAL dBt2 =  qBnd[5] - qInt[5];
    const REAL dphi =  qBnd[6] - qInt[6];
    const REAL dpsi =  qBnd[7] - qInt[7];

    // The four incoming amplitudes of that difference.
    const REAL w1 = dEt1 - c0*dBt2;
    const REAL w2 = dEt2 + c0*dBt1;
    const REAL w3 = dEn  - c0*dphi;
    const REAL w4 = dBn  - dpsi/c0;

    // Each contributes only to its own pair, through the incoming eigenvector:
    // (E_t1, B_t2) += w1 (1/2, -1/(2c)), (E_t2, B_t1) += w2 (1/2, +1/(2c)),
    // (E_n, phi)   += w3 (1/2, -1/(2c)), (B_n, psi)   += w4 (1/2, -c/2).
    const REAL cEn  = 0.5*w3,      cPhi = -0.5*w3/c0;
    const REAL cEt1 = 0.5*w1,      cBt2 = -0.5*w1/c0;
    const REAL cEt2 = 0.5*w2,      cBt1 =  0.5*w2/c0;
    const REAL cBn  = 0.5*w4,      cPsi = -0.5*w4*c0;

    // Back out of the face frame, and add to the interior state.
    qGhost[0] = qInt[0] + cEn*n0 - cEt1*n1;
    qGhost[1] = qInt[1] + cEn*n1 + cEt1*n0;
    qGhost[2] = qInt[2] + cEt2;
    qGhost[3] = qInt[3] + cBn*n0 - cBt1*n1;
    qGhost[4] = qInt[4] + cBn*n1 + cBt1*n0;
    qGhost[5] = qInt[5] + cBt2;
    qGhost[6] = qInt[6] + cPhi;
    qGhost[7] = qInt[7] + cPsi;
}

#endif // APMAXWELLCHARACTERISTICS_H
