#ifndef __wxadvectioneqn__
#define __wxadvectioneqn__

// WarpX includes
#include <wxcryptset.h>

// WarpX hyperbolic solver includes
#include <wxhyperboliceqn.h>

// std includes

template<typename REAL>
class WxAdvectionEqn : public WxHyperbolicEqn<REAL>
{
public:
    
/**
 * Construct new advection equation system
 */
    WxAdvectionEqn()
            : WxHyperbolicEqn<REAL>("advectionEqn") {
    }

/**
 * Number of equations in system
 */
    unsigned meqn() const {
        return 1;
    }

/**
 * Number of waves in Riemann solver
 */
    unsigned mwave() const {
        return 1;
    }

/**
 * Setup advection equation using supplied cryptset
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

/** Flux Jacobian **/
    void fluxJacobian(unsigned d, REAL *x, REAL *q, REAL *qaux, REAL **f);

/**
 * For an equation of the form Q_t = div(F) return the RHS
 *
 * @param N spatial order
 * @param geometry - [dr/dx, ds/dx, dr/dy, ds/dy]
 * @param normals - [nx1, nx2, nx3, ny1, ny2, ny3]
 * @param Dr - differentiation matrix in r-direction
 * @param Ds - differentiation matrix in s-direction
 * @param q [in] Conserved variable
 * @param dq [in] Conserved variable jump at element interface
 * @param rhs [out] div(F) evaluation
 */
   virtual void RHS(unsigned N, REAL *geometry, REAL *normals, WxpDGGeometry<REAL> *quad, REAL *q, REAL *dq, REAL *rhs);

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
 * Riemann solver for system where first order solution is used when
 * negative density/pressure is encountered
 *
 * @param d  direction along which Riemann solution is required
 * @param ql [in] Conserved variables at left edge
 * @param qr [in] Conserved variables at right edge
 * @param ql1 [in] Conserved variables at left edge using 1st order
 * @param qr1 [in] Conserved variables at right edge using 1st order
 * @param qauxl [in] Auxillary variables at left edge
 * @param qauxr [in] Auxillary variables at right edge
 * @param df [in] Jump across interface: this should be split into waves and fluctuations
 * @param wave [out] An meqn by mwaves array of waves in the system
 * @param s [out] Wave speeds
 * @param amdq [out] Negative (left) going fluctuations
 * @param apdq [out] Positive (right) going fluctuations
 */
    void
    rplimit(unsigned d, REAL *ql, REAL *qr, REAL *ql1, REAL *qr1, REAL *qauxl, REAL *qauxr, REAL *df, REAL **wave, REAL *s, REAL *amdq, REAL *apdq);

/**
 * Solve the transverse Riemann solve.
 *
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



    /**
     * Solve the transverse Riemann solve for the correction wave.
     *
     * @param td Direction in which the transverse correction split is desired
     * @param d Direction in which the data is oriented
     * @param ql [in] Conserved variables at left edge
     * @param qr [in] Conserved variables at right edge
     * @param s [in] Which correction wave is applied
     * @param bms Down going portion of correction wave
     * @param bps Up going portion of correction wave
     *
     *
     */

    void
    rptc(unsigned td, unsigned d, REAL *ql, REAL* qr,
    	REAL *s, REAL* bms, REAL* bps);


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
    REAL _ux, _uy, _uz;
};

#endif //  __wxadvectioneqn__
