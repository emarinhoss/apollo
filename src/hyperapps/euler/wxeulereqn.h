#ifndef __wxeulereqn__
#define __wxeulereqn__

// WarpX includes
#include <wxcryptset.h>
#include <wxmath.h>

// WarpX hyperbolic solver includes
#include <wxhyperboliceqn.h>

// std includes

template<typename REAL>
class WxEulerEqn : public WxHyperbolicEqn<REAL>
{
  public:
    
/**
 * Construct new Euler equation system
 */
    WxEulerEqn()
      : WxHyperbolicEqn<REAL>("eulerEqn") {
    }

/**
 * Number of equations in system
 */
    unsigned meqn() const {
      return 5;
    }

/**
 * Number of waves in Riemann solver
 */
    unsigned mwave() const {
      return 3;
    }

/**
 * Setup Euler equations using supplied cryptset
 */
    void setup(const WxCryptSet& wxc);

/**
 * Flux function
 *
 * @param d direction along which flux is required
 * @param x Coordinate at which flux is required
 * @param q [in] Conserved variable
 * @param qaux [in] Auxillary variable
 * @param f [out] Flux along direction 'd'
 */
    void flux(unsigned d, REAL *x, REAL *q, REAL *qaux, REAL *f);

/**
 * Flux Jacobian function
 *
 * @param d direction along which flux is required
 * @param x Coordinate at which flux is required
 * @param q [in] Conserved variable
 * @param qaux [in] Auxillary variable
 * @param f [out] Flux along direction 'd'
 */
    void fluxJacobian(unsigned d, REAL *x, REAL *q, REAL *qaux, REAL **f);

/**
 * Riemann solver for system
 *
 * @param td  direction along which transverse Riemann split is required
 * @param d  direction along which normal Riemann solution is required
 * @param ql [in] Conserved variables at left edge
 * @param qr [in] Conserved variables at right edge
 * @param qauxl [in] Auxillary variables at left edge
 * @param qauxr [in] Auxillary variables at right edge
 * @param df [in] Jump across interface: this should be split into waves and fluctuations
 * @param wave [out] An meqn by mwaves array of waves in the system
 * @param s [out] Wave speeds
 * @param amdq [out] Negative (left) going fluctuations
 * @param apdq [out] Positive (right) going fluctuations
 * @param soc [out] second order corrections
 */
    void
    rp(unsigned d, REAL *ql, REAL *qr, REAL *qauxl, REAL *qauxr, REAL *df, REAL **wave, REAL *s, REAL *amdq, REAL *apdq);
    
    void
    rpt(unsigned td, unsigned d, REAL *ql, REAL* qr, REAL *amdq, REAL* bmamdq, REAL* bpamdq, REAL *apdq, REAL* bmapdq, REAL* bpapdq);

    void
    rptc(unsigned td, unsigned d, REAL *ql, REAL* qr,REAL *soc, REAL* bms, REAL* bps);
/**
 * Compute eigensystem of flux Jacobian.
 *
 * @param d direction
 * @param q conserved variables
 * @param ev Eigenvalues
 * @param lev Left eigenvectors
 * @param rev Right eigenvectors
 */
    void eigenSystem(unsigned d, REAL *q, REAL *ev, REAL **lev, REAL **rev);

  private:
    REAL _gas_gamma, _minPres;
    bool _efix;
};

#endif //  __wxeulereqn__
