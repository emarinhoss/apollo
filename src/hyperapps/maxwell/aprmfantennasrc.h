#ifndef APRMFANTENNASRC_H
#define APRMFANTENNASRC_H

// WarpX includes
#include <wxlogger.h>
#include <wxlogstream.h>

// WarpX hyperbolic solver includes
#include <wxhyperbolicsrc.h>

// std includes
#include <cmath>
#include <string>

/**
 * A rotating-magnetic-field ANTENNA, as an axial current density in the domain.
 *
 * This is the drive the RMF literature actually uses: a pair of saddle coils
 * carrying axial currents with a cos(theta) distribution, rotated at omega by
 * driving the two coils in quadrature. The classical result is that such a
 * winding produces a UNIFORM transverse field inside itself - which is what
 * makes "rotating magnetic field" the right name for what the plasma sees -
 * and that field is an OUTCOME here, not a prescription.
 *
 * That distinction is the whole point of this class. Apollo's other RMF drives
 * (twoFluidSimplifiedRMFBC and friends) impose the vacuum field at the plasma
 * edge, which forbids the plasma from screening it: the antenna cannot be
 * loaded, the penetration threshold in B_omega is not the experiment's, and the
 * absorbed power cannot be computed. Driving from a current instead lets the
 * plasma's response reach the antenna and change the field there.
 * See docs/rmf-frc-model-assessment.md.
 *
 * THE CURRENT
 *
 *     J_z(r, theta, t) = K(t) f(r) cos(theta - omega t - phase)
 *
 * with f(r) a raised-cosine radial profile of width `width` centred on
 * `radius`, normalised so that the integral of f dr is one. A delta-function
 * sheet is not representable on a mesh, and a real antenna is a distributed
 * winding in any case; `width` should span at least two or three cells.
 *
 * The source it contributes is dE_z/dt += -J_z/epsilon0, matching the sign and
 * normalisation of the plasma current source (WxCurrentSrc, `Kind = currents`),
 * which computes -q n u / epsilon0. Only the z-component of E is touched, so
 * the deck's OutRange must be the single index of E_z.
 *
 * THE FIELD IT PRODUCES, which is how the amplitude is set
 *
 * Outside the plasma the transverse field is B = curl(A_z zhat) with
 * grad^2 A_z = -mu0 J_z. A single current shell of strength k at radius s,
 * inside a perfectly conducting cylinder at r = b, gives a uniform interior
 * field mu0 k (1 - s^2/b^2)/2 - the image current in the shell cancels part of
 * it. Superposing the shells of the profile,
 *
 *     B_uniform = (mu0/2) * integral of K f(s) (1 - s^2/b^2) ds
 *
 * and with b absent (free space, or a transparent outer boundary) the bracket
 * is one. This class evaluates that integral at setup and reports it. Rather
 * than make the user solve it backwards, the deck states the field it wants:
 *
 *     B_rmf              the uniform transverse field to produce [T]
 *     radius             antenna radius r_c [m]
 *     width              radial extent of the winding [m]
 *     conductor_radius   optional; radius of the perfectly conducting return
 *                        path. Omit for free space.
 *
 * and the amplitude K is solved for. `total_current` in the setup log is the
 * peak current per coil, integral of J_z over half the winding, which is the
 * number to compare against an experiment.
 *
 * ROTATION
 *
 * cos(theta - omega t - phase) rotates for any phase: the phase is an offset in
 * the rotation, never a change of polarisation. (The sibling boundary condition
 * applied its phase to one Cartesian component only, which silently turned the
 * drive into a linearly polarised oscillating field; writing the drive in this
 * form makes that class of mistake unrepresentable.)
 *
 * The envelope 1 - exp(-t/rise_time) ramps the current on, as the boundary
 * conditions do, so the run does not start with a step in J.
 */
template<class REAL>
class ApRMFAntennaSrc : public WxHyperbolicSrc<REAL>
{
  public:

    ApRMFAntennaSrc()
      : WxHyperbolicSrc<REAL>("maxwellRMFAntenna") {
    }

    void setup(const WxCryptSet& wxc) {
      WxHyperbolicSrc<REAL>::setup(wxc);

      _pi = 3.141592653589793;

      REAL freq  = wxc.template get<REAL>("frequency");
      _omega     = 2*_pi*freq;
      _phase     = wxc.template get<REAL>("phase");
      _rise      = wxc.template get<REAL>("rise_time");
      _rc        = wxc.template get<REAL>("radius");
      _width     = wxc.template get<REAL>("width");
      _epsilon0  = wxc.template get<REAL>("epsilon0");
      REAL mu0   = wxc.template get<REAL>("mu0");
      REAL Btgt  = wxc.template get<REAL>("B_rmf");

      // Free space unless the deck names a conducting return path.
      _b = wxc.has("conductor_radius")
             ? wxc.template get<REAL>("conductor_radius")
             : 0.0;

      _rmin = _rc - 0.5*_width;
      _rmax = _rc + 0.5*_width;

      // Reject geometries the shell superposition below does not describe,
      // rather than dividing by zero and driving the run with an infinity.
      // Each of these is a plausible deck typo.
      WxLogger *log = WxLogger::get("apollo-root.console");
      if (_width <= 0. || _rise <= 0. || _rc <= 0.) {
        log->error("*** maxwellRMFAntenna: radius, width and rise_time must all "
                   "be positive ***\n");
        exit(1);
      }
      if (_rmin <= 0.) {
        log->error("*** maxwellRMFAntenna: the winding reaches the axis "
                   "(radius <= width/2); it must be an annulus ***\n");
        exit(1);
      }
      if (_b > 0. && _b <= _rmax) {
        log->error("*** maxwellRMFAntenna: conductor_radius is inside the "
                   "winding; the antenna must sit inside its return path ***\n");
        exit(1);
      }

      // Shape integral: integral of f(s) (1 - s^2/b^2) ds over the winding.
      // Simpson on an odd number of points; the integrand is smooth, and the
      // profile's own normalisation is computed the same way so that the two
      // quadratures cancel to the extent that they are both inexact.
      const int N = 2001;
      REAL norm = 0., shape = 0.;
      REAL ds = (_rmax - _rmin)/(N-1);
      for (int i = 0; i < N; ++i) {
        REAL s = _rmin + i*ds;
        REAL w = (i == 0 || i == N-1) ? 1. : (i % 2 ? 4. : 2.);
        REAL fs = profile(s);
        REAL screen = (_b > 0.) ? (1. - s*s/(_b*_b)) : 1.;
        norm  += w*fs;
        shape += w*fs*screen;
      }
      norm  *= ds/3.;
      shape *= ds/3.;

      // profile() is not normalised on its own; _fnorm makes integral f ds = 1.
      _fnorm = 1./norm;
      REAL shapeNorm = shape/norm;

      // B_uniform = mu0/2 * K * shapeNorm  ->  solve for K.
      _K = 2.*Btgt/(mu0*shapeNorm);

      // -J/epsilon0, folded once so src() is a few multiplies.
      _KbyEps = _K/_epsilon0;

      // Peak current per coil: J_z integrated over the half of the winding
      // where the cosine is positive, at peak envelope. The theta integral of
      // cos over (-pi/2, pi/2) is 2, and the radial integral of r f(r) is
      // approximately r_c since f is normalised and narrow.
      _coilCurrent = 2.*_K*_rc;

      WxLogStream infoStrm = log->getInfoStream();
      infoStrm << "maxwellRMFAntenna: B_rmf = " << Btgt << " T requires K = "
               << _K << " A/m (screening factor " << shapeNorm
               << ", peak current per coil " << _coilCurrent << " A)"
               << std::endl;
    }

/**
 * Axial antenna current, as a source for dE_z/dt.
 *
 * @param n Length of q array (unused: the antenna does not depend on the state)
 * @param tx tx[0] is time and tx[1..3] are the spatial location
 * @param q Conserved variables at tx
 * @param qaux Auxiliary variables at tx
 * @param s Source term; s[0] is dE_z/dt
 */
    bool src(unsigned n, REAL *tx, REAL *q, REAL *qaux, REAL *s) {

      REAL t = tx[0];
      REAL x = tx[1];
      REAL y = tx[2];
      REAL r = sqrt(x*x + y*y);

      // Always write s[0]. WxHyperbolicSrc::compSource keeps _outValues as a
      // member and does sfull[idx] += _outValues[i] without clearing it, so a
      // src() that returns early leaves the PREVIOUS quadrature point's value
      // in place and adds it here. (ApRMFSrc in aprmfsrc.h has exactly that
      // defect.) An early return is never safe in this interface.
      s[0] = 0.0;

      if (r < _rmin || r > _rmax)
        return true;

      REAL envelope = 1. - exp(-t/_rise);
      // atan2, not atan(y/x): the latter folds the half-planes together and
      // would drive two half-antennas in phase instead of in opposition.
      REAL theta = atan2(y, x);

      // J_z/eps0, not J_z: _KbyEps folds the division in, since the only use
      // of the current here is as a source for dE_z/dt = ... - J_z/eps0.
      REAL JzByEps = _KbyEps*envelope*_fnorm*profile(r)
                     *cos(theta - _omega*t - _phase);

      s[0] = -JzByEps;

      return true;
    }

  private:

/**
 * Unnormalised raised cosine over the winding: smooth, and zero with zero
 * derivative at both edges, so the current does not step at the antenna's
 * boundary where the mesh is coarse.
 */
    REAL profile(REAL r) const {
      if (r < _rmin || r > _rmax) return 0.;
      return 0.5*(1. - cos(2.*_pi*(r - _rmin)/(_rmax - _rmin)));
    }

    REAL _omega, _phase, _rise, _rc, _width, _rmin, _rmax;
    REAL _b, _epsilon0, _K, _KbyEps, _fnorm, _coilCurrent, _pi;
};

#endif // APRMFANTENNASRC_H
