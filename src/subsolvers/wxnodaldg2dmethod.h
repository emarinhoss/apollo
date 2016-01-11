#ifndef WXNODALDG2DMETHOD_H
#define WXNODALDG2DMETHOD_H

// WarpX hyper includes
#include <wxhyperboliceqnset.h>
#include <wxhyperbolicsrcset.h>
#include <apareaintegralset.h>

// WarpX subsolver includes
#include <apsolver.h>
#include "apsubsolver.h"
#include <wxgridbc.h>
#include <wxnodaldglimiter.h>

// WarpX lib includes
#include <wxfunction.h>
#include <wxobject.h>
#include <wxnodaldggeometry2d.h>
#include <wxcubature2d.h>

// Petsc includes
#include <petscts.h>
#include <petscdmplex.h>
#include <petscsf.h>

// std includes
#include <map>
#include <string>
#include <vector>

template <typename REAL>
class WxNodalDG2dMethod : public ApSubSolver<REAL>
{
  public:
    WxNodalDG2dMethod()
            : ApSubSolver<REAL>("nodalDG2d") {
    }

/** Destructor */
    virtual ~WxNodalDG2dMethod();

/**
 * Setup subsolver object using supplied cryptset
 *
 * @param wxc Cryptset to use for setting
 */
    void setup(const WxCryptSet& wxc, DM dm);

/**
 * Initialize the subsolver: this is called after the setup() and
 * before the step() methods.
 */
    void init(PetscReal newDt, Vec out);

/**
 * Step the solver by given time step. If this step failed and if the
 * parent solver is running in variable time-stepping mode, it will be
 * called again with a smaller time-step. Otherwise, on failure the
 * simulation will be aborted.
 *
 * @param dt time step
 * @return false if step failed.
 */
    WxStepperStatus<REAL> step(REAL t, REAL dt, Vec in, Vec out);

/**
 * Return number of total equations in system
 *
 * @return Number of equations in system
 */
    unsigned totalEqns() const {
      return _eqnSet.totalEqns();
    }

/**
 * Returns the format of the dataStruture in this subsolver
 */
    std::vector<WxAny> getDataStructure(){
        return _dataStruct;
    }

/**
 *  Apply boundary conditions
 */
    void applyBc(int bcNum, REAL *xc, REAL *nx, REAL *q, REAL *qaux, REAL *AreaInts, REAL *qBC);

/**
 *  Apply boundary conditions
 */
    void applyLimiter(Vec Qin, Vec Qlimited);

 /**
 * Compute the RHS using CG spatial discretization.
 *
 * @param dt Time step to take
 * @param q Conserved variables
 * @param rhs Right hand side using
 */
   //WxStepperStatus<REAL> computeRhs(REAL dt, WxArray<REAL>& q, WxArray<REAL> &src);

/**
 * Checks to see if any the norm of a vector is infinity of not-a-number
 *
 * @param f - input petsc vector
 */

   PetscErrorCode isInfinityOrNAN(Vec f, std::string location);

private:

/**
 * Calculate the faces and cell geometry
 */

/** Time step to use */
  REAL _dt;
/** Array to modify the directions when doing cartesian v/s radial */
  unsigned _dirs[3];
/** Polynomial Order  */
  int _polyOrder;
/** cfl number */
  REAL _cfl;

/** Set of hyperbolic equations to solve */
  WxHyperbolicEqnSet<REAL> _eqnSet;
/** Set of source terms in equation system */
  WxHyperbolicSrcSet<REAL> _srcSet;
/** Set of area integral terms in equation system */
  ApAreaIntegralSet<REAL> _areaSet;

/** Arrays for passing to and fro from Reimann solver */
  // jumps, cons. var in left and right of edge i
  REAL *_qM, *_qP, *_fM, *_fP, *_gM, *_gP, *_qauxM, *_qauxP, *_numericalFLux;
  REAL *_apdq, *_amdq; // fluctuations
  REAL *_sx, *_sy; // speeds
  REAL **_wave; // waves
  REAL *_src; // source
  REAL *_areaInts, *_AgregateAreaIntegral; // area Integrals
/** Arrays for passing to and from from Reimann solver */
  REAL *_df;
/** Equations and waves */
  unsigned _meqn, _mwave;
/** element length, dx **/
  REAL _dx;
/** need to access status from step function */
  WxStepperStatus<REAL> _status;
  WxFunction<REAL>* _initFunc; // initial condition to use
  std::vector<WxAny> _dataStruct;
/** DG matrices, node coordinates and connectivities */
  wxNodalDGgeometry2D<REAL> *_geom;
/** Cubature information */
  WxCubature2d<REAL> *_cub;

/** Stuff need for the DG calculations */
  DM _dm;
/** local vectors with the data */
  Vec locU;

/** list of BC and limter subsolvers */
  std::vector<std::string> _bcSubSolvers;
  std::vector<std::string> _limiterSubSolvers;
  bool _haveLimiter;

/** Filtering variables **/
  REAL *_filterdiag, *_filterMatrix, _cutoff;
  int _orderSP;

};

#endif // WXNODALDG2DMETHOD_H
