/**
 * The rotating-magnetic-field boundary condition.
 *
 * Two properties that follow from the name of the thing:
 *
 *  1. The applied transverse field ROTATES. Its magnitude is constant over an
 *     RMF period and its direction sweeps through 2*pi. The `phase` parameter
 *     offsets where in that rotation the field starts; it must not change the
 *     polarisation. Before this test existed, `phase` was applied to only one
 *     of the two components, so phase = pi/2 - the value in the heavyIons
 *     deck - produced a linearly polarised field along a fixed axis whose
 *     magnitude swung between 0 and sqrt(2)*B_rmf. That is an oscillating
 *     field, a different experiment with its own boundary condition.
 *
 *  2. The axial electric field written at the boundary is the one Faraday's
 *     law induces from that rotating field. For E = E_z zhat and a spatially
 *     uniform B_perp(t), dB/dt = -curl E gives
 *
 *         dE_z/dx = +dB_y/dt ,     dE_z/dy = -dB_x/dt
 *
 *     so E_z = x dB_y/dt - y dB_x/dt, which rotates with the field. The
 *     previous expression carried the correct radial amplitude with the
 *     azimuthal dependence dropped and satisfied neither identity. E_z drives
 *     the oscillating axial currents that produce the azimuthal torque, so its
 *     angular structure is the torque's structure.
 *
 * Both are checked numerically against the shipped class, at the parameters of
 * examples/unstructuredDG/multifluid/rmf_frc.
 */

#include <wxcryptset.h>

#include "multifluid/aptwofluidsimplifiedrmfbc.h"

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

/** setup() is protected on the boundary-condition hierarchy. */
struct TestableRMFBC : public APTwoFluidSimplifiedRMFBC<double>
{
    void configure(const WxCryptSet& wxc) { this->setup(wxc, NULL); }
};

// Parameters of the hydrogen deck.
const double PLASMA_RADIUS = 0.030;
const double FLUX_RADIUS   = 0.035;
const double B_AXIAL       = 60.0e-4;
const double B_RMF         = 50.0e-4;
const double FREQ          = 7.95e5;
const double RISE          = 3.0e-8;

TestableRMFBC* makeBC(double phase)
{
    // %e, not %g: the deck parser types a literal without a decimal point as an
    // integer, and get<REAL> on it throws std::bad_cast with no context.
    // %.17g renders 7.95e5 as "795000", which is exactly that trap.
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
    // The cryptset must outlive setup(); the class copies what it needs.
    static WxCryptSet* held = NULL;
    delete held;
    held = new WxCryptSet(stream);
    TestableRMFBC* bc = new TestableRMFBC();
    bc->configure(held->getSet("coil"));
    return bc;
}

/** Evaluate the boundary state at (x, y, t). Returns E_z, B_x, B_y, B_z. */
void sample(TestableRMFBC& bc, double x, double y, double t, double out[4])
{
    double xc[3] = {t, x, y};
    double nx[2] = {1.0, 0.0};
    double q[18]  = {0};
    double qBC[18] = {0};
    double areaInts[18] = {0};
    // A tenuous but valid plasma state; the field part of the BC does not use
    // it, and the fluid part must not divide by zero.
    q[0] = 1.0e-12; q[5] = 1.0e-12; q[15] = B_AXIAL;
    bc.applyToArray(xc, nx, q, NULL, areaInts, qBC);
    out[0] = qBC[12]; out[1] = qBC[13]; out[2] = qBC[14]; out[3] = qBC[15];
}

/** The applied field must have constant magnitude and sweep a full turn. */
void testFieldRotates(double phase, const char* label)
{
    std::printf("\napplied transverse field, phase = %s\n", label);
    TestableRMFBC* bc = makeBC(phase);

    const double period = 1.0 / FREQ;
    // Sample well after the switch-on. The envelope 1 - exp(-t/rise) is still
    // climbing at t = 20*rise (its residual is exp(-20) = 2e-9, which shows up
    // as a spread in |B| of exactly that size); by 45*rise it is saturated to
    // below double precision, so any remaining variation in |B| is a genuine
    // departure from circular polarisation.
    const double t0 = 45.0 * RISE;
    const int N = 128;

    double mn = 1e300, mx = -1e300, sweep = 0.0, prev = 0.0;
    for (int k = 0; k <= N; ++k) {
        double f[4];
        sample(*bc, PLASMA_RADIUS, 0.0, t0 + k * period / N, f);
        const double mag = std::hypot(f[1], f[2]);
        mn = std::min(mn, mag); mx = std::max(mx, mag);
        const double ang = std::atan2(f[2], f[1]);
        if (k) {
            double d = ang - prev;
            while (d >  M_PI) d -= 2*M_PI;
            while (d < -M_PI) d += 2*M_PI;
            sweep += d;
        }
        prev = ang;
    }

    char detail[256];
    std::snprintf(detail, sizeof detail,
                  "|B| in [%.6e, %.6e], spread %.2e of the maximum",
                  mn, mx, (mx - mn) / mx);
    check((mx - mn) / mx < 1e-13, "magnitude is constant over an RMF period", detail);

    std::snprintf(detail, sizeof detail,
                  "swept %.4f rad in one period (want +-2*pi = %.4f)",
                  sweep, 2*M_PI);
    check(std::fabs(std::fabs(sweep) - 2*M_PI) < 1e-6,
          "direction sweeps exactly one full turn", detail);

    std::snprintf(detail, sizeof detail, "|B| = %.6e, B_rmf = %.6e", mx, B_RMF);
    check(std::fabs(mx - B_RMF) / B_RMF < 1e-6,
          "saturated magnitude equals B_rmf", detail);
    delete bc;
}

/** E_z must be the field Faraday's law induces from the applied B_perp. */
void testFaraday(double phase, const char* label)
{
    std::printf("\ninduced axial E-field, phase = %s\n", label);
    TestableRMFBC* bc = makeBC(phase);

    // Away from the axis and from the rise transient, so the central
    // differences below are well conditioned.
    const double x0 = 0.011, y0 = -0.007, t0 = 4.0e-7;
    const double h = 1.0e-6;      // metres
    const double dt = 1.0e-13;    // seconds

    double px[4], mx_[4], py[4], my[4], pt[4], mt[4];
    sample(*bc, x0 + h, y0, t0, px);
    sample(*bc, x0 - h, y0, t0, mx_);
    sample(*bc, x0, y0 + h, t0, py);
    sample(*bc, x0, y0 - h, t0, my);
    sample(*bc, x0, y0, t0 + dt, pt);
    sample(*bc, x0, y0, t0 - dt, mt);

    const double dEz_dx = (px[0] - mx_[0]) / (2*h);
    const double dEz_dy = (py[0] - my[0]) / (2*h);
    const double dBx_dt = (pt[1] - mt[1]) / (2*dt);
    const double dBy_dt = (pt[2] - mt[2]) / (2*dt);

    const double scale = std::max(std::fabs(dBx_dt), std::fabs(dBy_dt));
    char detail[256];

    std::snprintf(detail, sizeof detail,
                  "dEz/dx = %+.6e,  dBy/dt = %+.6e,  relative error %.2e",
                  dEz_dx, dBy_dt, std::fabs(dEz_dx - dBy_dt) / scale);
    check(std::fabs(dEz_dx - dBy_dt) / scale < 1e-6,
          "dEz/dx = +dBy/dt", detail);

    std::snprintf(detail, sizeof detail,
                  "dEz/dy = %+.6e, -dBx/dt = %+.6e,  relative error %.2e",
                  dEz_dy, -dBx_dt, std::fabs(dEz_dy + dBx_dt) / scale);
    check(std::fabs(dEz_dy + dBx_dt) / scale < 1e-6,
          "dEz/dy = -dBx/dt", detail);

    // The dropped-azimuth form was proportional to r with no angular
    // dependence, so it took the same value at every point on a circle.
    // The correct one does not.
    double a1[4], a2[4];
    const double rr = 0.02;
    sample(*bc, rr, 0.0, t0, a1);
    sample(*bc, 0.0, rr, t0, a2);
    std::snprintf(detail, sizeof detail,
                  "E_z at (r,0) = %+.6e, at (0,r) = %+.6e", a1[0], a2[0]);
    check(std::fabs(a1[0] - a2[0]) > 1e-9 * std::max(std::fabs(a1[0]), 1e-30),
          "E_z varies around the boundary circle", detail);
    delete bc;
}

/** The flux conserver must reproduce the bias field for a uniform column. */
void testFluxConserver()
{
    std::printf("\naxial flux conserver\n");
    TestableRMFBC* bc = makeBC(0.0);

    // A plasma still carrying the full bias flux must see the bias field at
    // the wall: pi a^2 B_axial in, (pi b^2 B_axial - pi a^2 B_axial)
    // /(pi(b^2-a^2)) = B_axial out.
    double xc[3] = {1.0e-9, PLASMA_RADIUS, 0.0};
    double nx[2] = {1.0, 0.0};
    double q[18] = {0}, qBC[18] = {0}, areaInts[18] = {0};
    q[0] = 1.0e-12; q[5] = 1.0e-12; q[15] = B_AXIAL;
    areaInts[15] = M_PI * PLASMA_RADIUS * PLASMA_RADIUS * B_AXIAL;
    bc->applyToArray(xc, nx, q, NULL, areaInts, qBC);

    char detail[256];
    std::snprintf(detail, sizeof detail, "B_z(a) = %.8e, B_axial = %.8e",
                  qBC[15], B_AXIAL);
    check(std::fabs(qBC[15] - B_AXIAL) / B_AXIAL < 1e-12,
          "unperturbed column sees the bias field at the wall", detail);

    // Expelling half the plasma flux must raise the wall field.
    areaInts[15] = 0.5 * M_PI * PLASMA_RADIUS * PLASMA_RADIUS * B_AXIAL;
    bc->applyToArray(xc, nx, q, NULL, areaInts, qBC);
    std::snprintf(detail, sizeof detail,
                  "half the flux expelled -> B_z(a) = %.6e (was %.6e)",
                  qBC[15], B_AXIAL);
    check(qBC[15] > B_AXIAL,
          "expelled flux is compressed into the annulus", detail);
    delete bc;
}

} // namespace

int main()
{
    std::printf("RMF boundary condition\n");
    std::printf("======================\n");

    testFieldRotates(0.0, "0");
    testFieldRotates(M_PI / 2.0, "pi/2 (the heavyIons deck)");
    testFaraday(0.0, "0");
    testFaraday(M_PI / 2.0, "pi/2");
    testFluxConserver();

    std::printf("\n%d checks, %d failures\n", checks, failures);
    return failures == 0 ? 0 : 1;
}
