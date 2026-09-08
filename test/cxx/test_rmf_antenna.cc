/**
 * The RMF antenna source term (maxwellRMFAntenna, ApRMFAntennaSrc).
 *
 * The antenna is the point of Phase 1 of docs/rmf-frc-model-assessment.md: the
 * rotating field is produced by a current rather than written onto the plasma
 * edge, so that the plasma can screen it and the antenna can be loaded. That
 * only means anything if the current the class emits really is the classical
 * RMF winding, so this checks it against the field theory says it produces,
 * without going near a mesh.
 *
 * The central check inverts the class. A z-directed surface current k cos(theta
 * - psi) on a shell of radius s, inside a perfectly conducting cylinder at
 * r = b, produces a UNIFORM transverse field inside itself,
 *
 *     A_z = (mu0 k/2)(1 - s^2/b^2) r cos(theta - psi)
 *     B   = (dA_z/dy, -dA_z/dx) = alpha (sin psi, -cos psi)
 *
 * the image current in the shell accounting for the (1 - s^2/b^2). Sampling
 * J_z from src() on a polar grid, extracting its m = 1 Fourier component shell
 * by shell and superposing gives the field the winding would make - which must
 * equal the B_rmf the deck asked for. That is a test of setup() and src()
 * together, along a path that shares no arithmetic with either.
 */

#include <wxcryptset.h>

#include "maxwell/aprmfantennasrc.h"

#include <cmath>
#include <cstdio>
#include <sstream>
#include <string>
#include <vector>

namespace {

int failures = 0;
int checks = 0;

void check(bool ok, const std::string& what, const std::string& detail = "")
{
    ++checks;
    std::printf("  %s  %s\n", ok ? "ok  " : "FAIL", what.c_str());
    if (!detail.empty())
        std::printf("        %s\n", detail.c_str());
    if (!ok) ++failures;
}

// Parameters of examples/unstructuredDG/multifluid/rmf_frc/antenna/frc2d.pin.
const double MU0     = M_PI*4.0e-7;
const double LIGHT   = 3.0e6;
const double EPS0    = 1.0/(MU0*LIGHT*LIGHT);
const double B_RMF   = 50.0e-4;
const double FREQ    = 7.95e5;
const double OMEGA   = 2.0*M_PI*FREQ;
const double RISE    = 3.0e-8;
const double COIL_R  = 0.040;
const double COIL_W  = 0.008;
const double WALL_R  = 0.060;

/** setup() is public on WxHyperbolicSrc, but OutRange etc. must be present. */
ApRMFAntennaSrc<double>* makeSrc(double phase, bool withConductor)
{
    // %e, not %g: the deck parser types a literal with no decimal point as an
    // integer, and get<REAL> on it throws std::bad_cast with no context.
    char deck[2048];
    std::snprintf(deck, sizeof deck,
        "<top>\n<ant>\n"
        "  Type = WxHyperbolicSrc\n"
        "  Kind = maxwellRMFAntenna\n"
        "  OutRange = [12]\n"
        "  frequency = %.17e\n"
        "  B_rmf = %.17e\n"
        "  phase = %.17e\n"
        "  rise_time = %.17e\n"
        "  radius = %.17e\n"
        "  width = %.17e\n"
        "  epsilon0 = %.17e\n"
        "  mu0 = %.17e\n"
        "%s"
        "</ant>\n</top>\n",
        FREQ, B_RMF, phase, RISE, COIL_R, COIL_W, EPS0, MU0,
        withConductor ? "  conductor_radius = 6.00000000000000000e-02\n" : "");
    std::istringstream stream(deck);
    static WxCryptSet* held = NULL;
    delete held;
    held = new WxCryptSet(stream);
    ApRMFAntennaSrc<double>* src = new ApRMFAntennaSrc<double>();
    src->setup(held->getSet("ant"));
    return src;
}

/** J_z at (x, y, t), recovered from the source the class emits. */
double currentAt(ApRMFAntennaSrc<double>& a, double x, double y, double t)
{
    double tx[3] = {t, x, y};
    double q[1] = {0.0};
    // Deliberately dirty, to catch a src() that returns without writing:
    // WxHyperbolicSrc::compSource keeps this array between calls and only ever
    // adds it to the state, so a stale value is silently added everywhere.
    double s[1] = {-12345.0};
    a.src(1, tx, q, NULL, s);
    // src() emits dE_z/dt = -J_z/eps0.
    return -s[0]*EPS0;
}

/**
 * The uniform interior field the sampled current would produce.
 *
 * Shell by shell: c(s) and d(s) are the cos and sin parts of the m = 1
 * component of J_z on the shell, so the shell carries k cos(theta - psi) with
 * (k cos psi, k sin psi) = (c, d), and contributes alpha (sin psi, -cos psi)
 * with alpha = (mu0/2)(1 - s^2/b^2) k ds. In components that is simply
 * (mu0/2)(1 - s^2/b^2) (d, -c) ds.
 */
void fieldFromCurrent(ApRMFAntennaSrc<double>& a, double t, bool withConductor,
                      double* bx, double* by)
{
    const int NR = 4001, NTH = 720;
    const double r0 = COIL_R - 0.75*COIL_W, r1 = COIL_R + 0.75*COIL_W;
    const double dr = (r1 - r0)/(NR - 1), dth = 2.0*M_PI/NTH;
    *bx = 0.0; *by = 0.0;
    for (int i = 0; i < NR; ++i) {
        const double s = r0 + i*dr;
        double c = 0.0, d = 0.0;
        for (int k = 0; k < NTH; ++k) {
            const double th = k*dth;
            const double J = currentAt(a, s*std::cos(th), s*std::sin(th), t);
            c += J*std::cos(th)*dth;
            d += J*std::sin(th)*dth;
        }
        c /= M_PI; d /= M_PI;   // m = 1 Fourier amplitudes
        const double screen = withConductor ? (1.0 - s*s/(WALL_R*WALL_R)) : 1.0;
        // Trapezoid in r: the integrand vanishes at both ends of this window.
        const double w = (i == 0 || i == NR-1) ? 0.5 : 1.0;
        *bx += 0.5*MU0*screen*d*dr*w;
        *by -= 0.5*MU0*screen*c*dr*w;
    }
}

/** The winding must reproduce the field the deck asked for. */
void testAmplitude(bool withConductor, const char* label)
{
    std::printf("\nfield produced by the winding, %s\n", label);
    ApRMFAntennaSrc<double>* a = makeSrc(0.0, withConductor);

    // Late enough that the switch-on envelope is saturated to well below the
    // quadrature error.
    const double t = 40.0*RISE;
    double bx, by;
    fieldFromCurrent(*a, t, withConductor, &bx, &by);
    const double mag = std::hypot(bx, by);

    char detail[256];
    std::snprintf(detail, sizeof detail, "|B| = %.6e, B_rmf = %.6e, relative error %.2e",
                  mag, B_RMF, std::fabs(mag - B_RMF)/B_RMF);
    check(std::fabs(mag - B_RMF)/B_RMF < 1e-4,
          "the current reproduces the requested B_rmf", detail);
    delete a;
}

/** Phase must offset the rotation, never change the polarisation. */
void testRotation()
{
    std::printf("\nrotation\n");
    ApRMFAntennaSrc<double>* a = makeSrc(0.0, true);

    const double t0 = 40.0*RISE, period = 1.0/FREQ;
    const int N = 24;
    double mn = 1e300, mx = -1e300, sweep = 0.0, prev = 0.0;
    for (int k = 0; k <= N; ++k) {
        double bx, by;
        fieldFromCurrent(*a, t0 + k*period/N, true, &bx, &by);
        const double m = std::hypot(bx, by);
        mn = std::min(mn, m); mx = std::max(mx, m);
        const double ang = std::atan2(by, bx);
        if (k) {
            double dd = ang - prev;
            while (dd >  M_PI) dd -= 2*M_PI;
            while (dd < -M_PI) dd += 2*M_PI;
            sweep += dd;
        }
        prev = ang;
    }

    char detail[256];
    std::snprintf(detail, sizeof detail, "|B| in [%.6e, %.6e], spread %.2e",
                  mn, mx, (mx - mn)/mx);
    check((mx - mn)/mx < 1e-6, "magnitude is constant over a period", detail);

    std::snprintf(detail, sizeof detail, "swept %.6f rad (want 2*pi = %.6f)",
                  std::fabs(sweep), 2*M_PI);
    check(std::fabs(std::fabs(sweep) - 2*M_PI) < 1e-6,
          "direction sweeps exactly one full turn", detail);
    delete a;
}

/** A phase offset must rotate the field by exactly that angle. */
void testPhaseIsAnOffset()
{
    std::printf("\nphase\n");
    const double t = 40.0*RISE;
    const double phase = M_PI/2.0;

    ApRMFAntennaSrc<double>* a0 = makeSrc(0.0, true);
    double bx0, by0;
    fieldFromCurrent(*a0, t, true, &bx0, &by0);
    delete a0;

    ApRMFAntennaSrc<double>* a1 = makeSrc(phase, true);
    double bx1, by1;
    fieldFromCurrent(*a1, t, true, &bx1, &by1);
    delete a1;

    const double m0 = std::hypot(bx0, by0), m1 = std::hypot(bx1, by1);
    char detail[256];
    std::snprintf(detail, sizeof detail, "|B| = %.6e with phase 0, %.6e with phase pi/2",
                  m0, m1);
    check(std::fabs(m0 - m1)/m0 < 1e-9,
          "phase leaves the magnitude alone (it is not a polarisation change)", detail);

    // Advancing the phase by p is the same as turning the clock back by p/omega,
    // so the field lags by exactly p.
    double d = std::atan2(by1, bx1) - std::atan2(by0, bx0);
    while (d >  M_PI) d -= 2*M_PI;
    while (d < -M_PI) d += 2*M_PI;
    std::snprintf(detail, sizeof detail, "rotated by %.9f rad, phase is %.9f", d, -phase);
    check(std::fabs(std::fabs(d) - phase) < 1e-9,
          "phase rotates the field by exactly that angle", detail);
}

/** Everything outside the winding must be written as an explicit zero. */
void testNoSourceOutsideTheWinding()
{
    std::printf("\ncurrent is confined to the winding\n");
    ApRMFAntennaSrc<double>* a = makeSrc(0.0, true);
    const double t = 40.0*RISE;

    // currentAt() seeds s[0] with -12345 before every call, so a src() that
    // returned early would show up here as that value rather than as zero.
    // This is the defect ApRMFSrc (maxwellRMFSrc) has: WxHyperbolicSrc keeps
    // _outValues between calls and adds it to the state unconditionally.
    double worst = 0.0;
    const double rs[] = {0.0, 0.005, 0.020, 0.030, 0.0355, 0.0445, 0.050, 0.059};
    for (unsigned i = 0; i < sizeof rs/sizeof rs[0]; ++i)
        for (int k = 0; k < 16; ++k) {
            const double th = 2.0*M_PI*k/16;
            worst = std::max(worst, std::fabs(
                currentAt(*a, rs[i]*std::cos(th), rs[i]*std::sin(th), t)));
        }

    char detail[256];
    std::snprintf(detail, sizeof detail, "largest |J_z| outside the winding: %.3e", worst);
    check(worst == 0.0, "J_z is exactly zero outside the winding", detail);

    // And non-zero inside it, so the check above is not passing vacuously.
    const double inside = std::fabs(currentAt(*a, COIL_R, 0.0, t));
    std::snprintf(detail, sizeof detail, "|J_z| at the winding centre: %.6e A/m^2", inside);
    check(inside > 1.0, "J_z is non-zero inside the winding", detail);
    delete a;
}

/** The current must be m = 1: equal and opposite across the axis. */
void testDipoleStructure()
{
    std::printf("\nm = 1 structure\n");
    ApRMFAntennaSrc<double>* a = makeSrc(0.0, true);
    const double t = 40.0*RISE;

    double worst = 0.0, scale = 0.0;
    for (int k = 0; k < 32; ++k) {
        const double th = 2.0*M_PI*k/32;
        const double p = currentAt(*a, COIL_R*std::cos(th), COIL_R*std::sin(th), t);
        const double m = currentAt(*a, -COIL_R*std::cos(th), -COIL_R*std::sin(th), t);
        worst = std::max(worst, std::fabs(p + m));
        scale = std::max(scale, std::fabs(p));
    }
    char detail[256];
    std::snprintf(detail, sizeof detail,
                  "max |J(theta) + J(theta+pi)| = %.3e against a peak of %.3e",
                  worst, scale);
    check(worst < 1e-9*scale, "J_z(theta+pi) = -J_z(theta)", detail);
    delete a;
}

/** The conducting shell must weaken the field by the image-current factor. */
void testConductorScreening()
{
    std::printf("\nconducting return path\n");
    const double t = 40.0*RISE;

    // Same requested B_rmf: the class must ask for MORE current when the wall
    // is there, by exactly 1/(1 - r_c^2/b^2) for a thin winding.
    ApRMFAntennaSrc<double>* free_ = makeSrc(0.0, false);
    ApRMFAntennaSrc<double>* walled = makeSrc(0.0, true);
    const double jFree = std::fabs(currentAt(*free_, COIL_R, 0.0, t));
    const double jWall = std::fabs(currentAt(*walled, COIL_R, 0.0, t));
    delete free_; delete walled;

    const double expected = 1.0/(1.0 - COIL_R*COIL_R/(WALL_R*WALL_R));
    const double got = jWall/jFree;
    char detail[256];
    std::snprintf(detail, sizeof detail,
                  "current ratio %.6f, thin-winding prediction %.6f", got, expected);
    // The winding is 8 mm wide, so the screening factor varies across it and
    // the ratio is only close to the thin-shell value.
    check(std::fabs(got - expected)/expected < 0.02,
          "the wall costs the antenna 1/(1 - r_c^2/b^2) in current", detail);
}

} // namespace

int main()
{
    std::printf("RMF antenna source\n");
    std::printf("==================\n");

    testAmplitude(true,  "inside a conducting wall");
    testAmplitude(false, "in free space");
    testRotation();
    testPhaseIsAnOffset();
    testNoSourceOutsideTheWinding();
    testDipoleStructure();
    testConductorScreening();

    std::printf("\n%d checks, %d failures\n", checks, failures);
    return failures == 0 ? 0 : 1;
}
