/**
 * Characteristic-consistent boundary states for perfectly-hyperbolic Maxwell.
 *
 * Phase 2 of docs/rmf-frc-model-assessment.md asked for boundary conditions
 * that impose only the incoming characteristics, on the grounds that writing a
 * whole prescribed field into the ghost state over-specifies a hyperbolic
 * system. This checks the machinery that does that - and also establishes the
 * result that changes what Phase 2 is for:
 *
 *   under the flux this code actually uses, with the cleaning speeds every deck
 *   actually sets, the over-specification makes no difference at all.
 *
 * WxPHMaxwellEqn::DGnumericalFlux is Lax-Friedrichs with
 * lambda = max(c, chi c, gamma c). At chi = gamma = 1 every eigenvalue of the
 * flux Jacobian has magnitude c, so |A| = c I and that flux is exactly the
 * upwind flux - which by construction takes the outgoing characteristics from
 * the interior and ignores whatever the ghost says about them. Algebraically
 * the flux depends on the ghost only through (A - c I) q_ghost, and (A - c I)
 * annihilates precisely the outgoing eigenvectors.
 *
 * So testAgainstOldGhost() below asserts the two constructions agree to
 * roundoff at chi = gamma = 1, and disagree by tens of percent as soon as
 * either is anything else. The correction matters for a deck that changes
 * DIVB_SPEED or DIVE_SPEED, and for the slope limiter, which consumes the
 * ghost state directly rather than through a Riemann solve.
 *
 * Everything here runs against the real WxPHMaxwellEqn - the Jacobian is
 * assembled by calling its flux() on unit vectors - so a change to the flux
 * that invalidated the decomposition would fail these rather than pass them.
 */

#include <wxcryptset.h>

#include "maxwell/wxphmaxwelleqn.h"
#include "maxwell/apmaxwellcharacteristics.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
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

const double C0 = 3.0e6;   // as the rmf_frc decks set it

WxPHMaxwellEqn<double>* makeEqn(double chi, double gamma)
{
    // %e, not %g: the deck parser types a literal with no decimal point as an
    // integer, and get<REAL> on it throws std::bad_cast with no context.
    char deck[512];
    std::snprintf(deck, sizeof deck,
        "<top>\n<maxwell>\n"
        "  Type = WxHyperbolicEqn\n  Kind = phMaxwellEqn\n"
        "  c0 = %.17e\n  gamma = %.17e\n  chi = %.17e\n"
        "</maxwell>\n</top>\n", C0, gamma, chi);
    std::istringstream stream(deck);
    static WxCryptSet* held = NULL;
    delete held;
    held = new WxCryptSet(stream);
    WxPHMaxwellEqn<double>* eqn = new WxPHMaxwellEqn<double>();
    eqn->setup(held->getSet("maxwell"));
    return eqn;
}

/** Flux Jacobian in direction d, assembled from the equation object itself. */
void jacobian(WxPHMaxwellEqn<double>& eqn, unsigned d, double A[8][8])
{
    for (int j = 0; j < 8; ++j) {
        double e[8] = {0,0,0,0,0,0,0,0}, f[8];
        e[j] = 1.0;
        eqn.flux(d, NULL, e, NULL, f);
        for (int i = 0; i < 8; ++i) A[i][j] = f[i];   // f is linear in q
    }
}

/** A deterministic spread of states and normals, so failures reproduce. */
double pseudo(int k)
{
    // A small LCG rather than rand(): the sequence must not depend on the
    // platform's libc for a failure to be reproducible from the output.
    static unsigned s = 12345u;
    if (k == 0) s = 12345u;
    s = 1664525u*s + 1013904223u;
    return 2.0*(double(s >> 8)/double(1u << 24)) - 1.0;
}

/**
 * The eight left eigenvectors the header derives, against the real Jacobian.
 */
void testEigenvectors()
{
    std::printf("\ncharacteristic variables against WxPHMaxwellEqn::flux\n");
    WxPHMaxwellEqn<double>* eqn = makeEqn(1.0, 1.0);
    double A[8][8];
    jacobian(*eqn, 0, A);   // face frame: normal along x

    // [E_x, E_y, E_z, B_x, B_y, B_z, phi, psi]
    struct { const char* name; double lam; double L[8]; } inv[] = {
      {"v1 = E_t1 + c B_t2", +C0, {0, 1, 0, 0, 0,  C0, 0, 0}},
      {"v2 = E_t2 - c B_t1", +C0, {0, 0, 1, 0, -C0, 0, 0, 0}},
      {"v3 = E_n + c phi",   +C0, {1, 0, 0, 0, 0, 0,  C0, 0}},
      {"v4 = B_n + psi/c",   +C0, {0, 0, 0, 1, 0, 0, 0, 1.0/C0}},
      {"w1 = E_t1 - c B_t2", -C0, {0, 1, 0, 0, 0, -C0, 0, 0}},
      {"w2 = E_t2 + c B_t1", -C0, {0, 0, 1, 0,  C0, 0, 0, 0}},
      {"w3 = E_n - c phi",   -C0, {1, 0, 0, 0, 0, 0, -C0, 0}},
      {"w4 = B_n - psi/c",   -C0, {0, 0, 0, 1, 0, 0, 0, -1.0/C0}},
    };

    double worst = 0.0;
    for (unsigned k = 0; k < sizeof inv/sizeof inv[0]; ++k) {
        double scale = 0.0, err = 0.0;
        for (int j = 0; j < 8; ++j) {
            double LA = 0.0;
            for (int i = 0; i < 8; ++i) LA += inv[k].L[i]*A[i][j];
            err = std::max(err, std::fabs(LA - inv[k].lam*inv[k].L[j]));
            scale = std::max(scale, std::fabs(inv[k].lam*inv[k].L[j]));
        }
        worst = std::max(worst, err/scale);
    }
    char detail[256];
    std::snprintf(detail, sizeof detail,
                  "worst |L A - lambda L| over the eight, relative: %.3e", worst);
    check(worst < 1e-14, "all eight are left eigenvectors of the real flux", detail);
    delete eqn;
}

/** At chi = gamma = 1 the Lax-Friedrichs flux is the exact upwind flux. */
void testLaxFriedrichsIsUpwind()
{
    std::printf("\nthe numerical flux at chi = gamma = 1\n");
    WxPHMaxwellEqn<double>* eqn = makeEqn(1.0, 1.0);
    double A[8][8];
    jacobian(*eqn, 0, A);

    // |A| = c I is what makes LF exact here. Check it as A*A = c^2 I, which
    // avoids an eigen-decomposition and is equivalent for a diagonalisable A
    // whose eigenvalues are all +-c.
    double worst = 0.0;
    for (int i = 0; i < 8; ++i)
        for (int j = 0; j < 8; ++j) {
            double s = 0.0;
            for (int k = 0; k < 8; ++k) s += A[i][k]*A[k][j];
            worst = std::max(worst, std::fabs(s - (i == j ? C0*C0 : 0.0)));
        }
    char detail[256];
    std::snprintf(detail, sizeof detail, "max |A^2 - c^2 I| = %.3e against c^2 = %.3e",
                  worst, C0*C0);
    check(worst < 1e-14*C0*C0,
          "A^2 = c^2 I, so |A| = c I and Lax-Friedrichs is upwind", detail);
    delete eqn;
}

/** Outgoing invariants must come from the interior, incoming from the boundary. */
void testGhostSplitsTheCharacteristics()
{
    std::printf("\nthe ghost state\n");
    double worstOut = 0.0, worstIn = 0.0;
    pseudo(0);
    for (int trial = 0; trial < 200; ++trial) {
        const double th = M_PI*pseudo(1);
        const double nx[2] = {std::cos(th), std::sin(th)};
        double qI[8], qB[8], qG[8];
        for (int i = 0; i < 8; ++i) { qI[i] = pseudo(1); qB[i] = pseudo(1); }
        maxwellCharacteristicGhost(nx, C0, qI, qB, qG);

        // Recompute the invariants of each state in the face frame.
        struct F { double En, Et1, Et2, Bn, Bt1, Bt2, phi, psi; };
        F f[3];
        const double* qs[3] = {qI, qB, qG};
        for (int s = 0; s < 3; ++s) {
            const double* q = qs[s];
            f[s].En  =  q[0]*nx[0] + q[1]*nx[1];
            f[s].Et1 = -q[0]*nx[1] + q[1]*nx[0];
            f[s].Et2 =  q[2];
            f[s].Bn  =  q[3]*nx[0] + q[4]*nx[1];
            f[s].Bt1 = -q[3]*nx[1] + q[4]*nx[0];
            f[s].Bt2 =  q[5];
            f[s].phi =  q[6];
            f[s].psi =  q[7];
        }
        // outgoing: ghost must match the interior
        worstOut = std::max(worstOut, std::fabs((f[2].Et1 + C0*f[2].Bt2) - (f[0].Et1 + C0*f[0].Bt2)));
        worstOut = std::max(worstOut, std::fabs((f[2].Et2 - C0*f[2].Bt1) - (f[0].Et2 - C0*f[0].Bt1)));
        worstOut = std::max(worstOut, std::fabs((f[2].En  + C0*f[2].phi) - (f[0].En  + C0*f[0].phi)));
        worstOut = std::max(worstOut, std::fabs((f[2].Bn  + f[2].psi/C0) - (f[0].Bn  + f[0].psi/C0)));
        // incoming: ghost must match the prescribed field
        worstIn = std::max(worstIn, std::fabs((f[2].Et1 - C0*f[2].Bt2) - (f[1].Et1 - C0*f[1].Bt2)));
        worstIn = std::max(worstIn, std::fabs((f[2].Et2 + C0*f[2].Bt1) - (f[1].Et2 + C0*f[1].Bt1)));
        worstIn = std::max(worstIn, std::fabs((f[2].En  - C0*f[2].phi) - (f[1].En  - C0*f[1].phi)));
        worstIn = std::max(worstIn, std::fabs((f[2].Bn  - f[2].psi/C0) - (f[1].Bn  - f[1].psi/C0)));
    }
    char detail[256];
    std::snprintf(detail, sizeof detail, "worst departure over 200 states and normals: %.3e", worstOut);
    check(worstOut < 1e-12*C0, "outgoing invariants are the interior's", detail);
    std::snprintf(detail, sizeof detail, "worst departure over 200 states and normals: %.3e", worstIn);
    check(worstIn < 1e-12*C0, "incoming invariants are the prescribed field's", detail);
}

/** Two limits the construction must satisfy exactly. */
void testLimits()
{
    std::printf("\nlimits\n");
    pseudo(0);
    double worstOpen = 0.0;
    for (int trial = 0; trial < 100; ++trial) {
        const double th = M_PI*pseudo(1);
        const double nx[2] = {std::cos(th), std::sin(th)};
        double q[8], qG[8];
        for (int i = 0; i < 8; ++i) q[i] = pseudo(1);
        maxwellCharacteristicGhost(nx, C0, q, q, qG);
        for (int i = 0; i < 8; ++i)
            worstOpen = std::max(worstOpen, std::fabs(qG[i] - q[i]));
    }
    char detail[256];
    std::snprintf(detail, sizeof detail, "worst |ghost - interior| = %.3e", worstOpen);
    check(worstOpen == 0.0,
          "prescribing the interior state degenerates to zero-gradient outflow", detail);

    // A field carrying only incoming characteristics, against a ZERO interior,
    // must pass through whole - the ghost's outgoing part is the interior's,
    // which is zero, and its incoming part is the prescribed one. (Against a
    // non-zero interior it would not: the ghost would carry that interior's
    // outgoing invariants, which is the entire point of the construction.)
    const double nx[2] = {1.0, 0.0};
    double qI[8] = {0,0,0,0,0,0,0,0}, qB[8], qG[8];
    // Purely incoming: E_t1 = -c B_t2 style pairs, i.e. v = 0 for all four.
    qB[0] = 1.0;  qB[6] = -qB[0]/C0;                 // E_n = -c phi  -> v3 = 0
    qB[3] = 2.0;  qB[7] = -qB[3]*C0;                 // B_n = -psi/c  -> v4 = 0
    qB[1] = 3.0;  qB[5] = -qB[1]/C0;                 // E_t1 = -c B_t2 -> v1 = 0
    qB[2] = 4.0;  qB[4] =  qB[2]/C0;                 // E_t2 = +c B_t1 -> v2 = 0
    maxwellCharacteristicGhost(nx, C0, qI, qB, qG);
    double worstIn = 0.0, scale = 0.0;
    for (int i = 0; i < 8; ++i) {
        worstIn = std::max(worstIn, std::fabs(qG[i] - qB[i]));
        scale = std::max(scale, std::fabs(qB[i]));
    }
    std::snprintf(detail, sizeof detail, "worst |ghost - prescribed| = %.3e against %.3e",
                  worstIn, scale);
    check(worstIn < 1e-9*scale,
          "a purely incoming prescribed field passes through unchanged", detail);
}

/**
 * The finding: at chi = gamma = 1 this changes nothing, elsewhere it changes
 * the boundary by tens of percent.
 */
void testAgainstOldGhost(double chi, double gamma, bool expectIdentical)
{
    std::printf("\nold ghost (the whole prescribed field) against the new one, "
                "chi = %g, gamma = %g\n", chi, gamma);
    WxPHMaxwellEqn<double>* eqn = makeEqn(chi, gamma);

    pseudo(0);
    double worst = 0.0, scale = 0.0;
    for (int trial = 0; trial < 300; ++trial) {
        const double th = M_PI*pseudo(1);
        double nrm[3] = {std::cos(th), std::sin(th), 1.0};
        double qI[8], qB[8], qG[8];
        for (int i = 0; i < 8; ++i) { qI[i] = pseudo(1); qB[i] = pseudo(1); }
        maxwellCharacteristicGhost(nrm, C0, qI, qB, qG);

        double fOld[8], fNew[8], sOld, sNew;
        eqn->DGnumericalFlux(nrm, qI, qB, fOld, &sOld);
        eqn->DGnumericalFlux(nrm, qI, qG, fNew, &sNew);
        for (int i = 0; i < 8; ++i) {
            worst = std::max(worst, std::fabs(fOld[i] - fNew[i]));
            scale = std::max(scale, std::fabs(fOld[i]));
        }
    }
    const double rel = worst/scale;
    char detail[256];
    std::snprintf(detail, sizeof detail,
                  "max |flux_old - flux_new| / max|flux| = %.3e", rel);
    if (expectIdentical)
        check(rel < 1e-13,
              "identical: Lax-Friedrichs is upwind here, so the ghost's "
              "outgoing part cannot matter", detail);
    else
        check(rel > 1e-2,
              "materially different: the over-specification is a real defect "
              "away from chi = gamma = 1", detail);
    delete eqn;
}

} // namespace

int main()
{
    std::printf("Maxwell characteristic boundary states\n");
    std::printf("======================================\n");

    testEigenvectors();
    testLaxFriedrichsIsUpwind();
    testGhostSplitsTheCharacteristics();
    testLimits();
    testAgainstOldGhost(1.0, 1.0, true);
    testAgainstOldGhost(1.5, 1.0, false);
    testAgainstOldGhost(1.0, 2.0, false);

    std::printf("\n%d checks, %d failures\n", checks, failures);
    return failures == 0 ? 0 : 1;
}
