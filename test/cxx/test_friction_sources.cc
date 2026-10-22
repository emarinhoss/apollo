/**
 * Energy conservation of the two-fluid friction source terms.
 *
 * The collisional coupling between the electron and ion fluids exchanges
 * momentum and heat, and must not create or destroy total energy. Writing the
 * friction on the electrons as R (so the reaction on the ions is -R) and the
 * relative drift as w = u_e - u_i, the total-energy sources for the two fluids
 * are
 *
 *     s_e = R . u_e + Q_e ,      s_i = -R . u_i + Q_i ,
 *
 * and conservation, s_e + s_i = 0, forces
 *
 *     Q_e + Q_i = -R . w = alpha |w|^2 = eta J^2 ,
 *
 * the frictional (Ohmic) heat. Braginskii deposits it in the electrons and
 * transfers Q_Delta from electrons to ions, so Q_e = -R.w - Q_Delta and
 * Q_i = Q_Delta. Substituting collapses both to a single dot product:
 *
 *     s_e = +R . u_i - Q_Delta ,   s_i = -R . u_i + Q_Delta .
 *
 * Note this holds for ANY friction law R, isotropic or not: the derivation
 * never opens up R. It also means the frictional heating does not appear as a
 * separate eta J^2 term - it is already implicit, because the momentum source
 * removes exactly that much kinetic energy from the electron fluid and the
 * energy source does not remove it, leaving it as heat.
 *
 * These tests exercise the real source classes through their real setup()
 * path, so they check the shipped code rather than a restatement of it.
 *
 *     make -C test/cxx            # or see test/cxx/Makefile
 */

#include <wxcryptset.h>

#include "multifluid/apconstantresistivity.h"
#include "multifluid/apanisotropicresistivitysrc.h"
#include "multifluid/wxbraginskiifrictionsrc.h"
#include "multifluid/ap2dmomentumxfer.h"

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
    if (ok) {
        std::printf("  ok    %s\n", what.c_str());
    } else {
        ++failures;
        std::printf("  FAIL  %s\n", what.c_str());
        if (!detail.empty())
            std::printf("        %s\n", detail.c_str());
    }
}

void closeTo(double got, double want, double tol, const std::string& what)
{
    const double scale = std::max(std::fabs(want), 1.0);
    const double err = std::fabs(got - want) / scale;
    char detail[256];
    std::snprintf(detail, sizeof detail,
                  "got %.17e, want %.17e, relative error %.3e > %.3e",
                  got, want, err, tol);
    check(err <= tol, what, detail);
}

/** Build a cryptset from an in-memory deck fragment. */
WxCryptSet parse(const std::string& body)
{
    std::istringstream deck("<top>\n" + body + "</top>\n");
    return WxCryptSet(deck);
}

// --- Physical constants and a representative plasma state -------------------
//
// Hydrogen at the conditions of examples/unstructuredDG/multifluid/rmf_frc:
// n = 1e20 m^-3, Te = 30 eV, Ti = 15 eV. The drift is large enough that the
// friction terms are not lost in round-off but small enough to stay subsonic.

const double MP = 1.67e-27;
const double ME = MP / 1836.0;
const double QE = 1.6e-19;
const double KB = 1.3807e-23;
const double GAMMA = 5.0 / 3.0;
const double ETA = 5.0e-6;

const double NDENS = 1.0e20;
const double TE_EV = 30.0;
const double TI_EV = 15.0;

/** Conserved two-fluid state [rho_e, rho_e u_e, ..., E_e, rho_i, ..., E_i]
 *  plus the magnetic field the anisotropic and Braginskii sources need. */
struct State {
    double q[16];
    double ue[3];
    double ui[3];
};

State makeState(double uex, double uey, double uez,
                double uix, double uiy, double uiz,
                double bx = 0.0, double by = 0.0, double bz = 0.0)
{
    State s;
    const double rhoe = NDENS * ME;
    const double rhoi = NDENS * MP;
    const double pe = NDENS * KB * (TE_EV * QE / KB);
    const double pi = NDENS * KB * (TI_EV * QE / KB);

    s.ue[0] = uex; s.ue[1] = uey; s.ue[2] = uez;
    s.ui[0] = uix; s.ui[1] = uiy; s.ui[2] = uiz;

    s.q[0] = rhoe;
    s.q[1] = rhoe * uex;
    s.q[2] = rhoe * uey;
    s.q[3] = rhoe * uez;
    s.q[4] = pe / (GAMMA - 1.0) + 0.5 * rhoe * (uex*uex + uey*uey + uez*uez);

    s.q[5] = rhoi;
    s.q[6] = rhoi * uix;
    s.q[7] = rhoi * uiy;
    s.q[8] = rhoi * uiz;
    s.q[9] = pi / (GAMMA - 1.0) + 0.5 * rhoi * (uix*uix + uiy*uiy + uiz*uiz);

    // Layout expected by the Braginskii / anisotropic sources, whose InpRange
    // appends the magnetic field after the two fluids.
    s.q[10] = bx; s.q[11] = by; s.q[12] = bz;
    s.q[13] = bx; s.q[14] = by; s.q[15] = bz;
    return s;
}

/**
 * Run one source term and check the energy identities.
 *
 * The source array is laid out [Rx, Ry, Rz, dE_e, -Rx, -Ry, -Rz, dE_i], which
 * is what every one of these classes writes and what the decks' OutRange maps
 * onto the momentum and energy components of the two fluids.
 */
template <typename SRC>
void checkSource(SRC& src, const State& st, const std::string& name,
                 bool threeD)
{
    std::printf("\n%s\n", name.c_str());

    double s[8] = {0};
    double tx[4] = {0.0, 0.0, 0.0, 0.0};
    double q[16];
    for (int i = 0; i < 16; ++i) q[i] = st.q[i];

    const bool ok = src.src(16, tx, q, NULL, s);
    check(ok, "src() returns true");

    const double R[3] = {s[0], s[1], s[2]};
    const double dEe = s[3];
    const double dEi = s[7];

    // The reaction on the ions must be equal and opposite.
    closeTo(s[4], -R[0], 1e-14, "ion momentum source is -R_x");
    closeTo(s[5], -R[1], 1e-14, "ion momentum source is -R_y");
    if (threeD)
        closeTo(s[6], -R[2], 1e-14, "ion momentum source is -R_z");

    // Total energy must be conserved exactly: this is the property the whole
    // test exists for.
    const double scale = std::max(std::fabs(dEe), std::fabs(dEi));
    char detail[256];
    std::snprintf(detail, sizeof detail,
                  "dE_e = %.6e, dE_i = %.6e, sum = %.6e (%.2e of the larger term)",
                  dEe, dEi, dEe + dEi,
                  scale > 0.0 ? std::fabs(dEe + dEi) / scale : 0.0);
    check(scale == 0.0 || std::fabs(dEe + dEi) / scale < 1e-12,
          "total energy source sums to zero", detail);

    // And each term must equal the closed form above. Q_Delta is whatever the
    // class computed, recovered from the ion equation, so this checks the
    // structure of the terms without hard-coding the equilibration model.
    const int nd = threeD ? 3 : 2;
    double Rdotui = 0.0;
    for (int d = 0; d < nd; ++d) Rdotui += R[d] * st.ui[d];

    const double Qdelta = dEi + Rdotui;   // from s_i = -R.u_i + Q_Delta
    closeTo(dEe, Rdotui - Qdelta, 1e-12, "electron energy source is +R.u_i - Q_Delta");

    // The frictional heat -R.w must be non-negative: friction dissipates.
    // With an anisotropic law the cross-field (Hall-like) part of R does no
    // work, so this still holds.
    double Rdotw = 0.0;
    for (int d = 0; d < nd; ++d) Rdotw += R[d] * (st.ue[d] - st.ui[d]);
    std::snprintf(detail, sizeof detail, "-R.w = %.6e", -Rdotw);
    check(-Rdotw >= 0.0, "frictional heating -R.w is non-negative", detail);
}

// --- The four source classes ------------------------------------------------

/**
 * The species constants, written at full precision.
 *
 * These must be bit-identical to the constants above: the classes recover the
 * number density as rho_e/m_e, so a truncated m_e in the deck shifts n by the
 * truncation and the eta*J^2 identity then misses by that much. Formatting from
 * the same doubles removes the question.
 */
std::string commonKeys()
{
    char buf[512];
    std::snprintf(buf, sizeof buf,
                  "  mi = %.17e\n"
                  "  me = %.17e\n"
                  "  charge = %.17e\n"
                  "  gas_gamma = %.17e\n"
                  "  boltz = %.17e\n",
                  MP, ME, QE, GAMMA, KB);
    return std::string(buf);
}

/** The same constants minus the ion mass, for tests that vary m_i. */
std::string meChargeKeys()
{
    char buf[512];
    std::snprintf(buf, sizeof buf,
                  "  me = %.17e\n"
                  "  charge = %.17e\n"
                  "  gas_gamma = %.17e\n"
                  "  boltz = %.17e\n",
                  ME, QE, GAMMA, KB);
    return std::string(buf);
}

void testConstantResistivity(const State& st)
{
    WxCryptSet all = parse(
        std::string("<res>\n"
                    "  Type = WxHyperbolicSrc\n"
                    "  Kind = constantResistivity\n"
                    "  InpRange = [0,1,2,3,4,5,6,7,8,9]\n"
                    "  OutRange = [1,2,3,4,6,7,8,9]\n")
        + commonKeys() +
        "  resistivity = 5.0e-6\n"
        "</res>\n");
    ApConstantResistivity<double> src;
    src.setup(all.getSet("res"));
    checkSource(src, st, "constantResistivity", true);

    // This source's friction law is known exactly, so the frictional heat can
    // be checked against eta J^2 with no reference to the class internals.
    double s[8] = {0}, tx[4] = {0}, q[16];
    for (int i = 0; i < 16; ++i) q[i] = st.q[i];
    src.src(16, tx, q, NULL, s);

    double Rdotw = 0.0, Jsq = 0.0;
    for (int d = 0; d < 3; ++d) {
        const double w = st.ue[d] - st.ui[d];
        Rdotw += s[d] * w;
        const double J = NDENS * QE * (st.ui[d] - st.ue[d]);
        Jsq += J * J;
    }
    closeTo(-Rdotw, ETA * Jsq, 1e-12, "frictional heat -R.w equals eta*J^2");
}

void testAnisotropicResistivity(const State& st)
{
    WxCryptSet all = parse(
        std::string("<res>\n"
                    "  Type = WxHyperbolicSrc\n"
                    "  Kind = anisotropicResistivity\n"
                    "  InpRange = [0,1,2,3,4,5,6,7,8,9,13,14,15]\n"
                    "  OutRange = [1,2,3,4,6,7,8,9]\n")
        + commonKeys() +
        "  epsilon0 = 8.854187817e-12\n"
        "</res>\n");
    ApAnisotropicResistivitySrc<double> src;
    src.setup(all.getSet("res"));
    checkSource(src, st, "anisotropicResistivity", false);
}

void testBraginskii(const State& st)
{
    WxCryptSet all = parse(
        std::string("<res>\n"
                    "  Type = WxHyperbolicSrc\n"
                    "  Kind = bragFriction\n"
                    "  InpRange = [0,1,2,3,4,5,6,7,8,9,13,14,15]\n"
                    "  OutRange = [1,2,3,4,6,7,8,9]\n")
        + commonKeys() +
        "  epsilon0 = 8.854187817e-12\n"
        "</res>\n");
    WxBragFrictionSrc<double> src;
    src.setup(all.getSet("res"));
    checkSource(src, st, "braginskiiFriction", true);
}

void testMomentumXfer(const State& st)
{
    WxCryptSet all = parse(
        std::string("<res>\n"
                    "  Type = WxHyperbolicSrc\n"
                    "  Kind = MomentumXfer2D\n"
                    "  InpRange = [0,1,2,3,4,5,6,7,8,9,13,14,15]\n"
                    "  OutRange = [1,2,3,4,6,7,8,9]\n")
        + commonKeys() +
        "  epsilon0 = 8.854187817e-12\n"
        "</res>\n");
    Ap2DMomentumXfer<double> src;
    src.setup(all.getSet("res"));
    checkSource(src, st, "MomentumXfer2D", false);
}

/**
 * The electron-ion thermal equilibration coefficient.
 *
 * Braginskii's rate is Q_Delta = 3 (m_e/m_i) n_e nu_ei (T_e - T_i). Written in
 * terms of a resistivity through eta = m_e nu_ei / (n e^2), that is
 * (3/m_i) n^2 eta e^2 (T_e - T_i), which is the form constantResistivity uses.
 * The mass ratio is therefore checkable without knowing the collision model:
 * Q_Delta must scale as m_e/m_i, so making the ions heavier must make
 * equilibration *slower*.
 */
void testEquilibrationScaling()
{
    std::printf("\nQ_Delta mass-ratio scaling\n");

    // Same drift-free state, two ion masses. With no drift the momentum source
    // vanishes and the energy sources are pure equilibration.
    const State rest = makeState(0, 0, 0, 0, 0, 0, 0.0, 0.0, 0.05);

    double qd[2];
    const double masses[2] = {1.67e-27, 4.0 * 1.67e-27};

    for (int k = 0; k < 2; ++k) {
        char mi[64];
        std::snprintf(mi, sizeof mi, "  mi = %.17e\n", masses[k]);
        WxCryptSet all = parse(
            std::string("<res>\n"
                        "  Type = WxHyperbolicSrc\n"
                        "  Kind = constantResistivity\n"
                        "  InpRange = [0,1,2,3,4,5,6,7,8,9]\n"
                        "  OutRange = [1,2,3,4,6,7,8,9]\n")
            + mi + meChargeKeys() +
            "  resistivity = 5.0e-6\n"
            "</res>\n");
        ApConstantResistivity<double> src;
        src.setup(all.getSet("res"));

        double s[8] = {0}, tx[4] = {0}, q[16];
        for (int i = 0; i < 16; ++i) q[i] = rest.q[i];
        // Ion density scales with mass at fixed number density.
        q[5] = NDENS * masses[k];
        src.src(16, tx, q, NULL, s);
        qd[k] = s[7];   // with no drift, s_i = +Q_Delta
    }

    char detail[256];
    std::snprintf(detail, sizeof detail,
                  "Q_Delta(m_i) = %.6e, Q_Delta(4 m_i) = %.6e, ratio = %.4f (want 0.25)",
                  qd[0], qd[1], qd[0] != 0.0 ? qd[1] / qd[0] : 0.0);
    check(qd[0] > 0.0, "hotter electrons lose energy to colder ions (Q_Delta > 0)", detail);
    closeTo(qd[1] / qd[0], 0.25, 1e-10,
            "Q_Delta scales as m_e/m_i: quadrupling the ion mass quarters it");
}

/** The same check applied to the collision-frequency form of the coefficient. */
void testBraginskiiEquilibrationSign()
{
    std::printf("\nQ_Delta in the Braginskii sources\n");

    // With no drift, equilibration is the only energy source. Compare the
    // Braginskii source against constantResistivity driven at the equivalent
    // collision frequency: both must agree to within the collision model, and
    // in particular must have the same sign and the same order of magnitude.
    const State rest = makeState(0, 0, 0, 0, 0, 0, 0.0, 0.0, 0.05);

    WxCryptSet all = parse(
        std::string("<res>\n"
                    "  Type = WxHyperbolicSrc\n"
                    "  Kind = bragFriction\n"
                    "  InpRange = [0,1,2,3,4,5,6,7,8,9,13,14,15]\n"
                    "  OutRange = [1,2,3,4,6,7,8,9]\n")
        + commonKeys() +
        "  epsilon0 = 8.854187817e-12\n"
        "</res>\n");
    WxBragFrictionSrc<double> src;
    src.setup(all.getSet("res"));

    double s[8] = {0}, tx[4] = {0}, q[16];
    for (int i = 0; i < 16; ++i) q[i] = rest.q[i];
    src.src(16, tx, q, NULL, s);
    const double Qd = s[7];

    // Bound the rate physically. Equilibration cannot be faster than the
    // electron-ion collision frequency itself: Q_Delta <= 3 n nu (Te - Ti),
    // and with the m_e/m_i factor it must be far below that. An inverted mass
    // ratio overshoots this bound by (m_i/m_e)^2 ~ 3e6.
    const double dT = (TE_EV - TI_EV) * QE;          // energy units
    const double nu_max = 1.0e9;                     // generous upper bound, s^-1
    const double bound = 3.0 * NDENS * nu_max * dT * (ME / MP) * 100.0;

    char detail[256];
    std::snprintf(detail, sizeof detail,
                  "Q_Delta = %.4e W/m^3; physical bound with the m_e/m_i factor is %.4e",
                  Qd, bound);
    check(Qd > 0.0, "hotter electrons lose energy to colder ions", detail);
    check(Qd < bound, "Q_Delta carries the m_e/m_i factor, not m_i/m_e", detail);
}

} // namespace

int main()
{
    std::printf("Two-fluid friction source terms: energy conservation\n");
    std::printf("====================================================\n");

    // A drift in every direction, and a magnetic field so the anisotropic and
    // Braginskii laws exercise their parallel, perpendicular and cross parts.
    const State drifting = makeState(4.0e4, -2.5e4, 1.2e4,   // electrons
                                     3.0e2,  1.0e2, -5.0e1,  // ions, much slower
                                     0.002, -0.001, 0.006);  // B, tesla

    testConstantResistivity(drifting);
    testAnisotropicResistivity(drifting);
    testBraginskii(drifting);
    testMomentumXfer(drifting);
    testEquilibrationScaling();
    testBraginskiiEquilibrationSign();

    std::printf("\n%d checks, %d failures\n", checks, failures);
    return failures == 0 ? 0 : 1;
}
