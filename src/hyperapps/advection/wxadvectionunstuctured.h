#ifndef WXADVECTIONUNSTUCTURED_H
#define WXADVECTIONUNSTUCTURED_H

// WarpX includes
#include <wxcryptset.h>

// WarpX hyperbolic solver includes
#include <wxhyperboliceqn.h>

// std includes

template<typename REAL>
class WxAdvectionUnstructuredEqn : public WxHyperbolicEqn<REAL>
{
  public:

/**
 * Construct new advection equation system
 */
    WxAdvectionUnstructuredEqn()
            : WxHyperbolicEqn<REAL>("advectionUnstEqn") {
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
 * Rotate data to local coordinate system defined by three unit
 * vectors. These vectors satisfy tan1 x tan2 = norm.
 *
 * @param norm Unit normal to surface.
 * @param tan1 Unit vector tangent to surface.
 * @param tan2 Unit vector tangent to surface. This is orthonormal to tan1.
 * @param vin Vector of data to rotate. This has meqn() components.
 * @param vout Vector of rotated data. This has meqn() components.
 */
    virtual void rotateToLocalFrame(REAL norm[3], REAL tan1[3], REAL tan2[3], REAL *vin, REAL *vout);

/**
 * Rotate data to global coordinate system. The local coordinate
 * system from which the data is to be rotated is defined by three
 * unit vectors. These vectors satisfy tan1 x tan2 = norm.
 *
 * @param norm Unit normal to surface.
 * @param tan1 Unit vector tangent to surface.
 * @param tan2 Unit vector tangent to surface. This is orthonormal to tan1.
 * @param vin Vector of data to rotate. This has meqn() components.
 * @param vout Vector of rotated data. This has meqn() components.
 */
    virtual void rotateToGlobalFrame(REAL norm[3], REAL tan1[3], REAL tan2[3], REAL *vin, REAL *vout);

/**
 * Rotate data to local coordinate system. The local coordinate system
 * lie entirely in the XY plain and hence only the normal to the
 * surface is provided.
 *
 * @param norm Unit normal to surface.
 * @param vin Vector of data to rotate. This has meqn() components.
 * @param vout Vector of rotated data. This has meqn() components.
 */
    virtual void rotateToLocalFrame(REAL norm[3], REAL *vin, REAL *vout);

/**
 * Rotate data to global coordinate system. The local coordinate
 * system from which the data is to be rotated is entirely in the XY
 * plain and hence only the normal to the surface is provided.
 *
 * @param norm Unit normal to surface.
 * @param vin Vector of data to rotate. This has meqn() components.
 * @param vout Vector of rotated data. This has meqn() components.
 */
    virtual void rotateToGlobalFrame(REAL norm[3], REAL *vin, REAL *vout);

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

#endif // WXADVECTIONUNSTUCTURED_H
