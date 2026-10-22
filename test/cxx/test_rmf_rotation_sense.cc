/**
 * Which way each RMF drive turns, and what that implies for field reversal.
 *
 * This is item 0 of Phase 3 in docs/rmf-frc-model-assessment.md, reduced to the
 * part that can be settled without compute. The assessment observes that the two
 * shipped RMF drives rotate in OPPOSITE senses while carrying the same +60 G
 * bias, so at most one of them is set up to form a field-reversed configuration.
 * A claim like that, left in prose, survives exactly until someone edits a sign.
 * Here it is executable, and it costs milliseconds rather than the tens of RMF
 * periods a simulation would need to show the same thing.
 *
 * THE SIGN CHAIN, stated in full because every step is a place to go wrong.
 *
 *   (x, y, z) is right-handed; +z is out of the r-theta plane, along the axis of
 *   the cylinder the 2-D run is a slice of. theta increases counter-clockwise
 *   from +x, and e_theta = (-sin theta, cos theta) points counter-clockwise.
 *
 *   1. A drive turning counter-clockwise drags electrons counter-clockwise. In
 *      the synchronous limit u_e,theta = zeta * omega * r with zeta -> 1; that
 *      limit is what the penetrated-state expression in the literature assumes.
 *   2. J = sum_s q_s n_s u_s, and the electron charge is negative, so
 *      counter-clockwise electrons carry a CLOCKWISE current: J_theta < 0.
 *   3. Inside an infinite cylinder each shell of azimuthal current contributes
 *      mu0 J_theta dr to the axial field, so
 *          B_z(0) = mu0 * integral of J_theta dr = -mu0 n e omega zeta a^2 / 2,
 *      which is NEGATIVE for zeta > 0.
 *   4. A negative driven B_z opposes a +z bias. That is field reversal.
 *
 *   So: forming an FRC against a +z bias requires a COUNTER-CLOCKWISE RMF. A
 *   clockwise one reinforces the bias instead, which is a different experiment.
 *
 * What is checked here is steps 1-4 as arithmetic, plus which sense each of the
 * two shipped drives actually has - measured from the classes themselves rather
 * than read off their source. The rotation sense is extracted the way a
 * diagnostic would: for the antenna, from the m = 1 Fourier component of the
 * current it emits; for the boundary condition, from the direction of the
 * transverse field it writes. Neither shares arithmetic with the class.
 *
 * The same sign rule is asserted from the Python side, on the diagnostic that
 * post-processes a run, in test/test_rmf_diagnostics.py. Both must agree.
 */

#include <wxcryptset.h>

#include "maxwell/aprmfantennasrc.h"
#include "multifluid/aptwofluidsimplifiedrmfbc.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <sstream>
#include <string>

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

// Shared parameters of the two rmf_frc decks.
const double MU0     = M_PI*4.0e-7;
const double LIGHT   = 3.0e6;
const double EPS0    = 1.0/(MU0*LIGHT*LIGHT);
const double B_RMF   = 50.0e-4;
const double B_AXIAL = 60.0e-4;         // the bias, +z
const double FREQ    = 7.95e5;
const double OMEGA   = 2.0*M_PI*FREQ;
const double RISE    = 3.0e-8;
const double COIL_R  = 0.036;
const double COIL_W  = 0.005;
const double PLASMA_RADIUS = 0.030;
const double FLUX_RADIUS   = 0.035;

const double Q_E   = 1.6e-19;
const double N_E   = 1.0e20;

/** Wrap into (-pi, pi], so that successive angles can be differenced. */
double wrap(double d)
{
    while (d >  M_PI) d -= 2.0*M_PI;
    while (d <= -M_PI) d += 2.0*M_PI;
    return d;
}

// ---------------------------------------------------------------------------
// The antenna source.
// ---------------------------------------------------------------------------

ApRMFAntennaSrc<double>* makeAntenna(double phase)
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
        "</ant>\n</top>\n",
        FREQ, B_RMF, phase, RISE, COIL_R, COIL_W, EPS0, MU0);
    std::istringstream stream(deck);
    static WxCryptSet* held = NULL;
    delete held;
    held = new WxCryptSet(stream);
    ApRMFAntennaSrc<double>* src = new ApRMFAntennaSrc<double>();
    src->setup(held->getSet("ant"));
    return src;
}

/**
 * The angle of the antenna's m = 1 current pattern at time t.
 *
 * J_z on the coil radius is k cos(theta - psi) for some k and psi; the Fourier
 * projections c = integral J cos theta and d = integral J sin theta give
 * (k cos psi, k sin psi), so psi = atan2(d, c). Taking it from the emitted
 * source rather than from the class's own expression is the point: this would
 * still measure the truth if the class computed the current some other way.
 */
double antennaPatternAngle(ApRMFAntennaSrc<double>& a, double t)
{
    const int NTH = 720;
    const double dth = 2.0*M_PI/NTH;
    double c = 0.0, d = 0.0;
    for (int k = 0; k < NTH; ++k) {
        const double th = k*dth;
        double tx[3] = {t, COIL_R*std::cos(th), COIL_R*std::sin(th)};
        double q[1] = {0.0};
        // Deliberately dirty: WxHyperbolicSrc::compSource keeps this array
        // between calls and only ever adds it, so a src() that returns without
        // writing silently reuses the previous point's value.
        double s[1] = {-12345.0};
        a.src(1, tx, q, NULL, s);
        const double J = -s[0]*EPS0;      // src emits dE_z/dt = -J_z/eps0
        c += J*std::cos(th)*dth;
        d += J*std::sin(th)*dth;
    }
    return std::atan2(d, c);
}

// ---------------------------------------------------------------------------
// The edge-driven boundary condition.
// ---------------------------------------------------------------------------

struct TestableRMFBC : public APTwoFluidSimplifiedRMFBC<double>
{
    void configure(const WxCryptSet& wxc) { this->setup(wxc, NULL); }
};

TestableRMFBC* makeBC(double phase)
{
    char deck[2048];
    std::snprintf(deck, sizeof deck,
        "<top>\n<coil>\n"
        "  Type = ApSubSolver\n"
        "  Kind = twoFluidSimplifiedRMFBC\n"
        "  frequency = %.17e\n"
        "  B_axial = %.17e\n"
        "  B_rmf = %.17e\n"
        "  phase = %.17e\n"
        "  rise_time = %.17e\n"
        "  plasma_radius = %.17e\n"
        "  flux_conserver_radius = %.17e\n"
        "</coil>\n</top>\n",
        FREQ, B_AXIAL, B_RMF, phase, RISE, PLASMA_RADIUS, FLUX_RADIUS);
    std::istringstream stream(deck);
    static WxCryptSet* held = NULL;
    delete held;
    held = new WxCryptSet(stream);
    TestableRMFBC* bc = new TestableRMFBC();
    bc->configure(held->getSet("coil"));
    return bc;
}

/** The direction of the transverse field the boundary condition writes. */
double bcFieldAngle(TestableRMFBC& bc, double t)
{
    double xc[3] = {t, PLASMA_RADIUS, 0.0};
    double nx[2] = {1.0, 0.0};
    double q[18] = {0}, qBC[18] = {0}, areaInts[18] = {0};
    q[0] = 1.0e-12; q[5] = 1.0e-12; q[15] = B_AXIAL;
    bc.applyToArray(xc, nx, q, NULL, areaInts, qBC);
    return std::atan2(qBC[14], qBC[13]);   // atan2(B_y, B_x)
}

/**
 * Total angle swept between t0 and t0 + n*dt, unwrapped.
 *
 * Signed, and that is the whole point: taking the magnitude would let a
 * handedness flip through, which is exactly the defect being guarded against.
 * dt must be well under half a period or the unwrapping is ambiguous.
 *
 * THE WINDOW MUST NOT BE A HALF-INTEGER NUMBER OF PERIODS. A standing
 * (linearly polarised) field has a direction angle that is fixed modulo pi and
 * jumps by pi each time its amplitude passes through zero - twice per period.
 * Over exactly half a period it therefore accumulates exactly +/-pi, which is
 * bit-identical to what a genuine rotation through half a turn produces. This
 * test measured over 0.5 of a period until someone checked. It now measures
 * 0.74, where a standing field can only give 0 or +/-pi and a rotation gives
 * 4.65 rad, and sweepIsMonotone below closes the case a coincidence could still
 * reach.
 */
template<class F>
double sweep(F angleAt, double t0, double dt, int n)
{
    double total = 0.0, prev = angleAt(t0);
    for (int i = 1; i <= n; ++i) {
        const double now = angleAt(t0 + i*dt);
        total += wrap(now - prev);
        prev = now;
    }
    return total;
}

/**
 * True if every step turns the same way by about the same amount.
 *
 * This is what separates rotating from standing without relying on the window
 * length at all. A rotation advances by omega*dt every step; a standing field
 * sits still and then jumps by pi. Requiring each increment to be within a
 * factor of two of the mean, and of the same sign, admits the first and refuses
 * the second whatever the window.
 */
template<class F>
bool sweepIsMonotone(F angleAt, double t0, double dt, int n)
{
    const double mean = sweep(angleAt, t0, dt, n)/n;
    if (std::fabs(mean) < 1e-12) return false;
    double prev = angleAt(t0);
    for (int i = 1; i <= n; ++i) {
        const double now = angleAt(t0 + i*dt);
        const double step = wrap(now - prev);
        prev = now;
        if (step*mean <= 0.0) return false;                 // turned back
        if (std::fabs(step) > 2.0*std::fabs(mean)) return false;   // jumped
        if (std::fabs(step) < 0.5*std::fabs(mean)) return false;   // stalled
    }
    return true;
}

// ---------------------------------------------------------------------------
// Step 3 of the sign chain, as arithmetic.
// ---------------------------------------------------------------------------

/**
 * Driven B_z on axis for synchronous electrons rotating with sense `sense`.
 *
 * sense = +1 counter-clockwise, -1 clockwise. Deliberately written out from
 * J_theta = -e n u_theta rather than from the -mu0 n e omega a^2/2 shorthand,
 * so that the electron charge sign appears explicitly where it can be checked.
 */
double drivenBzOnAxis(double sense, double zeta, double a)
{
    const double u_theta_coeff = sense*zeta*OMEGA;   // u_theta = coeff * r
    const double j_theta_coeff = -Q_E*N_E*u_theta_coeff;
    // B_z(0) = mu0 * integral from 0 to a of J_theta(r) dr
    return MU0*j_theta_coeff*a*a/2.0;
}

// ---------------------------------------------------------------------------

void testAntennaTurnsCounterClockwise(double phase, const char* label)
{
    std::printf("\nantenna current pattern, phase = %s\n", label);
    ApRMFAntennaSrc<double>* a = makeAntenna(phase);

    // Well past the ramp, so the envelope is not changing the amplitude much;
    // it cannot change the angle at all, but starting at t = 0 would sample a
    // vanishing current.
    const double t0 = 40.0*RISE, dt = 0.02/FREQ;
    const int n = 37;   // 0.74 of a period; see sweep()
    const double swept = sweep([&](double t) { return antennaPatternAngle(*a, t); },
                               t0, dt, n);
    const double expected = OMEGA*dt*n;

    char detail[256];
    std::snprintf(detail, sizeof detail,
                  "swept %+.6f rad in %.3e s, +omega*dt = %+.6f",
                  swept, dt*n, expected);
    check(swept > 0.0,
          "the antenna's current pattern turns COUNTER-CLOCKWISE", detail);
    check(std::fabs(swept - expected)/std::fabs(expected) < 1e-6,
          "and it turns at exactly omega", detail);
    check(sweepIsMonotone([&](double t) { return antennaPatternAngle(*a, t); },
                          t0, dt, n),
          "and it ROTATES rather than standing: every step turns the same way "
          "by the same amount", detail);
    delete a;
}

void testBoundaryTurnsClockwise(double phase, const char* label)
{
    std::printf("\nedge-driven boundary field, phase = %s\n", label);
    TestableRMFBC* bc = makeBC(phase);

    const double t0 = 40.0*RISE, dt = 0.02/FREQ;
    const int n = 37;   // 0.74 of a period; see sweep()
    const double swept = sweep([&](double t) { return bcFieldAngle(*bc, t); },
                               t0, dt, n);
    const double expected = -OMEGA*dt*n;

    char detail[256];
    std::snprintf(detail, sizeof detail,
                  "swept %+.6f rad in %.3e s, -omega*dt = %+.6f",
                  swept, dt*n, expected);
    check(swept < 0.0,
          "twoFluidSimplifiedRMFBC's field turns CLOCKWISE", detail);
    check(std::fabs(swept - expected)/std::fabs(expected) < 1e-6,
          "and it turns at exactly omega", detail);
    check(sweepIsMonotone([&](double t) { return bcFieldAngle(*bc, t); },
                          t0, dt, n),
          "and it ROTATES rather than standing: every step turns the same way "
          "by the same amount", detail);

    // A rotating field has constant magnitude; a standing one swings through
    // zero. Checked here as well as in test_rmf_boundary.cc because this file
    // is the one that claims to establish the drive is a ROTATING field.
    double lo = 1e300, hi = 0.0;
    for (int i = 0; i <= n; ++i) {
        double xc[3] = {t0 + i*dt, PLASMA_RADIUS, 0.0};
        double nx[2] = {1.0, 0.0};
        double q[18] = {0}, qBC[18] = {0}, areaInts[18] = {0};
        q[0] = 1.0e-12; q[5] = 1.0e-12; q[15] = B_AXIAL;
        bc->applyToArray(xc, nx, q, NULL, areaInts, qBC);
        const double mag = std::hypot(qBC[13], qBC[14]);
        lo = std::min(lo, mag); hi = std::max(hi, mag);
    }
    std::snprintf(detail, sizeof detail,
                  "|B| ranges over %.6e to %.6e T", lo, hi);
    check((hi - lo)/hi < 1e-6,
          "and its magnitude is constant, so it is circular and not linear "
          "polarisation", detail);
    delete bc;
}

void testTheTwoDrivesDisagree()
{
    std::printf("\nthe two drives against each other\n");
    ApRMFAntennaSrc<double>* a = makeAntenna(0.0);
    TestableRMFBC* bc = makeBC(0.0);

    const double t0 = 40.0*RISE, dt = 0.02/FREQ;
    const int n = 37;   // 0.74 of a period; see sweep()
    const double sa = sweep([&](double t) { return antennaPatternAngle(*a, t); },
                            t0, dt, n);
    const double sb = sweep([&](double t) { return bcFieldAngle(*bc, t); },
                            t0, dt, n);

    char detail[256];
    std::snprintf(detail, sizeof detail,
                  "antenna %+.6f rad, boundary condition %+.6f rad", sa, sb);
    check(sa*sb < 0.0,
          "they turn in OPPOSITE senses, so at most one can form an FRC "
          "against a given bias", detail);
    check(std::fabs(sa + sb) < 1e-9*std::fabs(sa),
          "and at the same rate, so this is handedness and nothing else",
          detail);
    delete a;
    delete bc;
}

void testSignChain()
{
    std::printf("\nrotation sense to driven axial field\n");

    const double bz_ccw = drivenBzOnAxis(+1.0, 1.0, PLASMA_RADIUS);
    const double bz_cw  = drivenBzOnAxis(-1.0, 1.0, PLASMA_RADIUS);

    char detail[256];
    std::snprintf(detail, sizeof detail,
                  "counter-clockwise: B_z = %+.4e T; clockwise: %+.4e T; "
                  "bias = %+.4e T", bz_ccw, bz_cw, B_AXIAL);
    check(bz_ccw < 0.0,
          "counter-clockwise electrons drive B_z < 0, opposing a +z bias",
          detail);
    check(bz_cw > 0.0,
          "clockwise electrons drive B_z > 0, reinforcing a +z bias", detail);

    // The magnitude, against the assessment's table: mu0 n e omega a^2/2 is
    // quoted there as 450 G for these parameters.
    const double gauss = std::fabs(bz_ccw)*1.0e4;
    std::snprintf(detail, sizeof detail,
                  "|B_z| at zeta = 1 is %.1f G; assessment section 3.2 says ~450 G, "
                  "%.1fx the %.0f G bias", gauss, gauss/(B_AXIAL*1.0e4), B_AXIAL*1.0e4);
    check(std::fabs(gauss - 452.0) < 5.0,
          "and its magnitude is the synchronous-limit value in the assessment",
          detail);
}

void testTheVerdict()
{
    std::printf("\nwhich deck is set up to form an FRC\n");

    ApRMFAntennaSrc<double>* a = makeAntenna(0.0);
    TestableRMFBC* bc = makeBC(0.0);
    const double t0 = 40.0*RISE, dt = 0.02/FREQ;
    const int n = 37;   // 0.74 of a period; see sweep()

    // The sense each drive turns, as +1 or -1, measured not assumed.
    const double sense_a = sweep([&](double t) { return antennaPatternAngle(*a, t); },
                                 t0, dt, n) > 0.0 ? +1.0 : -1.0;
    const double sense_b = sweep([&](double t) { return bcFieldAngle(*bc, t); },
                                 t0, dt, n) > 0.0 ? +1.0 : -1.0;

    // Push each through the sign chain, with a bias of +B_AXIAL.
    const double bz_a = drivenBzOnAxis(sense_a, 1.0, PLASMA_RADIUS);
    const double bz_b = drivenBzOnAxis(sense_b, 1.0, PLASMA_RADIUS);

    char detail[320];
    std::snprintf(detail, sizeof detail,
                  "antenna drives %+.4e T against a %+.4e T bias; "
                  "the boundary condition drives %+.4e T", bz_a, B_AXIAL, bz_b);
    check(bz_a*B_AXIAL < 0.0,
          "the ANTENNA deck's drive opposes its bias: set up for field reversal",
          detail);
    check(bz_b*B_AXIAL > 0.0,
          "the EDGE-DRIVEN deck's drive REINFORCES its bias: not field reversal",
          detail);
    delete a;
    delete bc;
}

}  // namespace

int main()
{
    std::printf("RMF rotation sense and the field-reversal sign chain\n");

    // Phase must offset the rotation, never change its sense. Both values are
    // shipped: 0 in the hydrogen decks, pi/2 in heavyIons.
    testAntennaTurnsCounterClockwise(0.0, "0");
    testAntennaTurnsCounterClockwise(M_PI/2.0, "pi/2");
    testBoundaryTurnsClockwise(0.0, "0");
    testBoundaryTurnsClockwise(M_PI/2.0, "pi/2");

    testTheTwoDrivesDisagree();
    testSignChain();
    testTheVerdict();

    std::printf("\n%d checks, %d failures\n", checks, failures);
    return failures ? 1 : 0;
}
