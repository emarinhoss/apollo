#ifndef __wxhyperboliceqn__
#define __wxhyperboliceqn__

// WarpX lib includes
#include <wxobject.h>
#include <wxcryptset.h>
#include <wxexcept.h>
#include <wxpdggeometry.h>

// std includes
#include <string>

/**
 * WxHyperbolicEqn represents a system of hyperbolic conservation laws.
 *
 * kappa*dq/dt + df/dx + dg/dy + dh/dz = s
 *
 * where kappa(x,y,z) is a capacity function, f, g, h are fluxes in
 * the X, Y and Z directions respectively and s is the sources.
 *
 */
template<typename REAL>
class WxHyperbolicEqn : public WxObject
{
  public:

/**
 * Construct new hyperbolic equation with supplied name.
 *
 * @param name name of equation system
 */
    WxHyperbolicEqn(const std::string& name);


    virtual ~WxHyperbolicEqn();

/**
 * Returns name of the equation system.
 *
 * @return Name of the equations system
 */
    std::string name() const;

/**
 * Returns number of equations in system.
 *
 * @return number of equations in system.
 */
    virtual unsigned meqn() const = 0;

/**
 * Returns number of variables to calculate gradients for equations in system.
 *
 * @return number of gradients for the equations in system.
 */
    virtual unsigned mgrads() const = 0;

/**
 * Returns number of waves in Riemann solver.
 *
 * @return Number of waves in Riemann solver.
 */
    virtual unsigned mwave() const = 0;

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
 * Sets spatial location of left cell center for use in Riemann problem
 * solver.
 *
 * @param ndim Number of dimensions.
 * @param x Location of left cell center.
 */
    void setCoordLeftCell(unsigned ndim, REAL *x);

/**
 * Sets spatial location of right cell center for use in Riemann
 * problem solver.
 *
 * @param ndim Number of dimensions.
 * @param x Location of right cell center.
 */
    void setCoordRightCell(unsigned ndim, REAL *x);

/**
 * Gets spatial location of left cell center for use in Riemann
 * problem solver.
 *
 * @param x Coordinates of left cell center.
 */
    void coordLeftCell(REAL *x);

/**
 * Gets spatial location of right cell center for use in Riemann
 * problem solver.
 *
 * @param x Coordinates of right cell center.
 */
    void coordRightCell(REAL *x);

/**
 * Flux function
 *
 * @param d direction along which flux is required
 * @param x Coordinate at which flux is required
 * @param q [in] Conserved variable
 * @param qaux [in] Auxillary variable
 * @param f [out] Flux along direction 'd'
 */
    virtual void flux(unsigned d, REAL *x, REAL *q, REAL *qaux, REAL *f);

/**
 * DG numerical Flux function
 *
 * @param normals [in] outward facing normal to a given face/edge
 * @param qM [in] conserved variables on the interior of the face
 * @param qP [in] conserved variables on the exterior of the face
 * @param f [out] Numerical Flux
 * @param maxSpeed [out] speed of fastest propagating wave
 */
    virtual void DGnumericalFlux(REAL *normals, REAL *qM, REAL *qP, REAL *nflux, REAL maxSpeed);

/**
 * DG limiter trigger variable
 *
 * @param qIn  [in] get all the variables for set equations system
 * @param qOut [out] returns the variable that is the trigger
 */
     virtual void DGLimiterTrigger(REAL *qIn, REAL *qOut);

/**
 * Flux Jacobian function
 *
 * @param d direction along which flux is required
 * @param x Coordinate at which flux is required
 * @param q [in] Conserved variable
 * @param qaux [in] Auxillary variable
 * @param f [out] Flux Jacobian along direction 'd'
 */
    virtual void fluxJacobian(unsigned d, REAL *x, REAL *q, REAL *qaux, REAL **f);

/**
 * Edge Flux function for DG gen geom
 *
 * @param d direction along which flux is required
 * @param x Coordinate at which flux is required
 * @param q [in] Conserved variable
 * @param qaux [in] Auxillary variable
 * @param f [out] Flux along direction 'd'
 */
    virtual void edgefluxgengeom(unsigned d, REAL *x, REAL *q, REAL *qaux, REAL *f);

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
   virtual void RHS(unsigned N, REAL *geometry, REAL *normals,
                            WxpDGGeometry<REAL> *quad, REAL *q, REAL *dq, REAL *rhs);

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
    virtual 
    void
    rp(unsigned d, REAL *ql, REAL *qr, 
      REAL *qauxl, REAL *qauxr, REAL *df, 
      REAL **wave, REAL *s, REAL *amdq, REAL *apdq) = 0;

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
    virtual
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

    virtual
    void
    rptc(unsigned td, unsigned d, REAL *ql, REAL* qr,
      REAL *s, REAL* bms, REAL* bps);

/**
 * Riemann solver for system with a built-in limiter
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
    virtual 
    void
    rplimit(unsigned d, REAL *ql, REAL *qr, REAL *ql1, REAL *qr1,
      REAL *qauxl, REAL *qauxr, REAL *df, 
      REAL **wave, REAL *s, REAL *amdq, REAL *apdq);

 /**
 * Compute eigensystem of flux Jacobian.
 *
 * @param d direction
 * @param q conserved variables
 * @param ev Eigenvalues
 * @param lev Left eigenvectors
 * @param rev Right eigenvectors
 */
    virtual 
    void
    eigenSystem(unsigned d, REAL *q, REAL *ev, REAL **lev, REAL **rev);


 /**
  * Set the normal, and tangential vectors
  * in a set of variables for access in a
  * local equation system method.  The set
  * is accomplished at the overall subsubsolver
  * level.
  * @param norm is the normal vector
  * @param tan1 is the first tangential vector
  * @param tan2 is the second tangential vector
  */

    void setFaceVectors(REAL norm[3], REAL tan1[3], REAL tan2[3]);

    /**
     * Get the normal, and tangential vectors
     * from a set of variables by a local equation
     * system method.
     * @param norm is the normal vector
     * @param tan1 is the first tangential vector
     * @param tan2 is the second tangential vector
     */

    void getFaceVectors(REAL norm[3], REAL tan1[3], REAL tan2[3]);

    /**
     * Set the indices
     * in a set of variables for access in a
     * local equation system method.  The set
     * is accomplished at the overall subsubsolver
     * level.
     * @param i is the dir[0] index
     * @param j is the dir[1] index
     * @param k is the dir[2] index
     */


    void setIndices(int idx[3]);

    /**
     * Set the indices
     * in a set of variables for access in a
     * local equation system method.  The set
     * is accomplished at the overall subsubsolver
     * level.
     * @param i is the dir[0] index
     * @param j is the dir[1] index
     * @param k is the dir[2] index
     */

    void getIndices(int idx[3]);

    /**
     * Calculate the primitive variables from the conserved variables
     *
     * @param qCons [in] conserved variables
     * @param qPrim [out] primitive variables
     *
     */

    virtual void primitiveVariables(REAL *qCons, REAL *qPrim);

    /**
     * Limit quantities using the Tu and Aliabadi limiter
     *
     * @param avgCons [in] average of the conservative variables
     * @param avgPrim [in] average of the primitive variables
     * @param dGrads  [in] gradients
     * @param limitedValues [out] positivity enforced values
     *
     */
    virtual void limiterTuAndAliabadi(REAL *avgCons, REAL *avgPrim, REAL *dGrads, REAL *limitedValues);


  private:
/** Name of equation system */
    std::string _name;
/** Left cell's cell-center coordinate */
    REAL _xl[3];
/** Right cell's cell-center coordinate */
    REAL _xr[3];
/** cell face vectors */
	REAL _norm[3], _tan1[3], _tan2[3];
/** cell indices 	 */
	int _idx[3];
};

#endif //  __wxhyperboliceqn__
