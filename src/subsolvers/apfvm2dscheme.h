#ifndef APFVM2DSCHEME_H
#define APFVM2DSCHEME_H

// WarpX hyper includes
#include <wxhyperboliceqnset.h>
#include <wxhyperbolicsrcset.h>

// WarpX subsolver includes
#include <apsolver.h>
#include "apsubsolver.h"

// WarpX lib includes
#include <wxfunction.h>
#include <wxobject.h>

// Petsc includes
#include <petscfv.h>
#include <petscts.h>

// std includes
#include <map>
#include <string>
#include <vector>

template <typename REAL>
class ApFVM2Dscheme : public ApSubSolver<REAL>
{
  public:
    ApFVM2Dscheme()
            : ApSubSolver<REAL>("fv2dscheme") {
    }

///** Constructor */
//    ApFVM2Dscheme();

/** Destructor */
    virtual ~ApFVM2Dscheme();

/**
 * Setup subsolver object using supplied cryptset
 *
 * @param wxc Cryptset to use for setting
 */
    void setup(const WxCryptSet& wxc);

/**
 * Initialize the subsolver: this is called after the setup() and
 * before the step() methods.
 */
    void init(Vec out);

/**
 * Step the solver by given time step. If this step failed and if the
 * parent solver is running in variable time-stepping mode, it will be
 * called again with a smaller time-step. Otherwise, on failure the
 * simulation will be aborted.
 *
 * @param dt time step
 * @return false if step failed.
 */
    WxStepperStatus<REAL> step(REAL dt, Vec in, Vec out);

/**
 * Return number of total equations in system
 *
 * @return Number of equations in system
 */
    unsigned totalEqns() const {
      return _eqnSet.totalEqns();
    }

 /**
 * Compute the RHS using CG spatial discretization.
 *
 * @param dt Time step to take
 * @param q Conserved variables
 * @param rhs Right hand side using
 */
   //WxStepperStatus<REAL> computeRhs(REAL dt, WxArray<REAL>& q, WxArray<REAL> &src);

private:

/** Time step to use */
  REAL _dt;
/** Array to modify the directions when doing cartesian v/s radial */
  unsigned _dirs[3];
/** Set of hyperbolic equations to solve */
  WxHyperbolicEqnSet<REAL> _eqnSet;
/** Set of source terms in equation system */
  WxHyperbolicSrcSet<REAL> _srcSet;
/** Arrays for passing to and fro from Reimann solver */
  // jumps, cons. var in left and right of edge i
  REAL *_ql, *_qr, *_fl, *_fr, *_qauxl, *_qauxr;
  REAL *_src; // source
/** Arrays for passing to and from from Reimann solver */
  REAL *_df;
  REAL *_s; // speed
/** Equations and waves */
  unsigned _meqn, _mwave;
/** element length, dx **/
  REAL _dx;
/** need to access status from step function */
  WxStepperStatus<REAL> _status;
/** data management and user-defined context */
  DM _dm;
  UserContext _usr;
  //PetscViewer viewer;
  Vec qvars;
  WxFunction<REAL>* _initFunc; // initial condition to use
  std::vector<std::string> _initArrays; // name of arrays to initialize

};

#endif // APFVM2DSCHEME_H
