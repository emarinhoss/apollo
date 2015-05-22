#ifndef __wxhyperboliceqnset__
#define __wxhyperboliceqnset__

// WarpX hyper includes
#include "wxhyperboliceqn.h"

// WarpX lib includes
#include <wxobject.h>
#include <wxexcept.h>

template<typename REAL>
class WxHyperbolicEqnSet : public WxObject
{
  public:
/** Create a new object */
    WxHyperbolicEqnSet() {}

/** Dtor */
    virtual ~WxHyperbolicEqnSet();

/**
 * Set with dirs vector which should adjust for radial problems. In
 * most cases dirs[] = {0,1,2}
 *
 * @param dirs[] Vector to adjust for radial problems.
 */    
    void setDirs(unsigned dirs[]);

/**
 * Set number of spatial dimensions for problem
 *
 * @param dim Spatial dimension
 */    
    void setDim(unsigned dim);

/**
 * Setup object using supplied crypset. This method should read in
 * parameters needed to initialize it.
 *
 * @param wxc Cryptset using which the object is set up.
 */
    virtual void setup(const WxCryptSet& wxc);

/**
 * Get total number to waves in equation set
 *
 * @return Number of waves in equation set
 */
    unsigned totalWaves() const {
      return _mwave;
    }

/**
 * Get total number of equations in set
 *
 * @return Number of equations in set
 */
    unsigned totalEqns() const {
      return _meqn;
    }

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
 * Computes flux from individual equations and assembles them into the
 * full flux for the equation system being solved
 *
 * @param d direction along which flux is required
 * @param x Spatial coordinates at which flux is required
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
 * @param nflux [out] Numerical Flux
 * @param maxSpeed [out] speed of fastest propagating wave
 */
    virtual void DGnumericalFlux(REAL *normals, REAL *qM, REAL *qP, REAL *nflux, REAL maxSpeed);

/**
 * Computes flux jacobian from individual equations and assembles them into the
 * full flux jacobian for the equation system being solved
 *
 * @param d direction along which flux is required
 * @param x Spatial coordinates at which flux is required
 * @param q [in] Conserved variable
 * @param qaux [in] Auxillary variable
 * @param f [out] Flux Jacobian along direction 'd'
 */
    void fluxJacobian(unsigned d, REAL *x, REAL *q, REAL *qaux, REAL **f);

/**
 * For general geometry DG implementation, only the normal component
 * of the flux are used. Therefore, this is direction independent. 
 * However, for the volume flux, the direction dependent flux function
 * is still needed
 *
 * @param d direction along which flux is required
 * @param x Spatial coordinates at which flux is required
 * @param q [in] Conserved variable
 * @param qaux [in] Auxillary variable
 * @param f [out] Flux along direction 'd'
 */
    void edgefluxgengeom(unsigned d, REAL *x, REAL *q, REAL *qaux, REAL *f);

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
    void calulateRHS(unsigned N, REAL *geometry, REAL *normals, WxpDGGeometry<REAL> *quad, REAL *q, REAL *dq, REAL *rhs);

/**
 * Solves the Reimann problem for individual systems and assembles
 * them into a larger system.
 */
    void
    riemann(unsigned d, REAL *xl, REAL *xr, REAL *ql, REAL *qr, REAL *qauxl, REAL *qauxr,
      REAL *df, REAL **wave, REAL *s, REAL *amdq, REAL *apdq);

/**
 * Solves the transvers Reimann problem for individual systems and
 * assembles them into a larger system.
 */
    void 
    riemannt(unsigned td, unsigned d, REAL *xl, REAL *xr, REAL *ql, REAL *qr,
      REAL *amdq, REAL* bmamdq, REAL* bpamdq,
      REAL *apdq, REAL* bmapdq, REAL* bpapdq);

/**
 * Solves the transvers Reimann problem for individual systems and
 * assembles them into a larger system.
 */
    void
    riemanntc(unsigned td, unsigned d, REAL *xl, REAL *xr, REAL *ql, REAL *qr,
      REAL *s, REAL* bms, REAL* bps);

/**
 * Solves the Reimann problem where negative densities and pressures
 * are replaced by first order solution
 */
    void
    riemannlimit(unsigned d, REAL *xl, REAL *xr, REAL *ql, REAL *qr, REAL *ql1, REAL *qr1, 
            REAL *qauxl, REAL *qauxr, REAL *df, REAL **wave, REAL *s, REAL *amdq, REAL *apdq);

/**
 * Computes the eigensystem for individual systems and assembles
 * them into a larger system.
 */
    void
    eigenSystem(unsigned d, REAL *q, REAL *ev, REAL **lev, REAL **rev);

    void
    setFaceVectors(REAL norm[3], REAL tan1[3], REAL tan2[3]);

    void
    setIndices(int idx[3]);

  private:

/** Private to prevent copying */
    WxHyperbolicEqnSet(const WxHyperbolicEqnSet<REAL>&);

/** Private to prevent assignment */
    WxHyperbolicEqnSet<REAL>& operator=(const WxHyperbolicEqnSet<REAL>&);

/** Total number of waves and equations */
    unsigned _mwave, _meqn;
/** List of equations to solve */
    std::vector<WxHyperbolicEqn<REAL>* > _eqnSys;
/** Direction vector to adjust for radial problems */
    unsigned _dirs[23];
/** Temporary space for waves for use in Riemann solver */
    REAL **_temp_wave, **_temp_flux;
    REAL **_temp_rev, **_temp_lev; // temporary space for right, left eigenvectors
/** Set number of dimensions */
    unsigned _ndim;
};

#endif // __wxhyperboliceqnset__
