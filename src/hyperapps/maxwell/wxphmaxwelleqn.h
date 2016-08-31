#ifndef __wxphmaxwelleqn__
#define __wxphmaxwelleqn__

// WarpX includes
#include <wxcryptset.h>
#include <wxmath.h>

// WarpX hyperbolic solver includes
#include <wxhyperboliceqn.h>

// std includes

template<typename REAL>
class WxPHMaxwellEqn : public WxHyperbolicEqn<REAL>
{
  public:
    
/**
 * Construct new PHMaxwell equation system
 */
    WxPHMaxwellEqn()
      : WxHyperbolicEqn<REAL>("phMaxwellEqn") {
    }

/**
 * Number of equations in system
 */
    unsigned meqn() const {
      return 8;
    }

/**
 * Number of waves in Riemann solver
 */
    unsigned mwave() const {
      return 6;
    }

/**
 * Setup PHMaxwell equations using supplied cryptset
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
 * DG numerical Flux function
 *
 * @param normals [in] outward facing normal to a given face/edge
 * @param qM [in] conserved variables on the interior of the face
 * @param qP [in] conserved variables on the exterior of the face
 * @param f [out] Numerical Flux
 * @param maxSpeed [out] speed of fastest propagating wave
 */
    void DGnumericalFlux(REAL *normals, REAL *qM, REAL *qP, REAL *nflux, REAL maxSpeed);

/**
 * Riemann solver for system
 *
 * @param d  direction along which Riemann solution is required
 * @param ql [in] Conserved variables at left edge
 * @param qr [in] Conserved variables at right edge
 * @param qauxl [in] Auxillary variables at left edge
 * @param qauxr [in] Auxillary variables at right edge
 * @param df [in] Jump across interface: this should be split into waves and fluctuations
 * @param wave [out] An meqn by mwaves array of waves in the system
 * @param s [out] Wave speeds
 * @param amdq [out] Negative (left) going fluctuations
 * @param apdq [out] Positive (right) going fluctuations
 */
    void
    rp(unsigned d, REAL *ql, REAL *qr, REAL *qauxl, REAL *qauxr, REAL *df, REAL **wave, REAL *s, REAL *amdq, REAL *apdq);

/**
 * Solve the transverse Riemann solve.
 *
 * @param td Direction in which the transverse split is desired
 * @param d Direction in which the data is oriented
 * @param ql [in] Conserved variables at left edge
 * @param qr [in] Conserved variables at right edge
 * @param imp [in] Which flux difference is supplied: imp=1 for A^-\Delta q and imp=2 for A^+ \Delta q
 * @param asdq [in] Flux differences in direction imp
 * @param bmasdq Down going portion of flux difference
 * @param bpasdq Up going portion of flux difference
 * 
 * For example, for a linear system q_t + Aq_x + Bq_y = 0, 
 * asdq = A^+ dq or A^- dq and this is then split into 
 * bmasdq = B^- asdq and bpasdq = B^+ asdq

 */
    void
    rpt(unsigned td, unsigned d, REAL *ql, REAL* qr,
      REAL *amdq, REAL* bmamdq, REAL* bpamdq,
      REAL *apdq, REAL* bmapdq, REAL* bpapdq);

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


/** Compute the primitive variables
 *
 * @param qCons [in] Conserved Variables
 * @param qPrim [out] Primitive Variables
 */
    void primitiveVariables(REAL *qCons, REAL *qPrim);


/**
 * Limit quantities using the Tu and Aliabadi limiter
 *
 * @param avgCons [in] average of the conservative variables
 * @param avgPrim [in] average of the primitive variables
 * @param dGrads  [in] gradients
 * @param limitedValues [out] positivity enforced values
 *
 */
    void limiterTuAndAliabadi(REAL *avgCons, REAL *avgPrim, REAL *dGrads, REAL *limitedValues);

  private:
    REAL _c0, _chi, _gamma;
};

#endif //  __wxphmaxwelleqn__
