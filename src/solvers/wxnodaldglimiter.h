#ifndef WXNODALDGLIMITER_H
#define WXNODALDGLIMITER_H

// WarpX subsolver includes
#include <apsubsolver.h>
#include "apsolver.h"
#include <wxnodaldggeometry2d.h>
#include <petscdmplex.h>

/**
 * Base class which applies boundary conditions to fields living on a
 * grid box.
 */
template <typename REAL>
class WxNodalDGLimiter : public ApSubSolver<REAL>
{
  public:

/**
 * Construct a new grid-bc object
 */
    WxNodalDGLimiter(const std::string& name)
      : ApSubSolver<REAL>(name) {
    }

    virtual ~WxNodalDGLimiter() {
    }

/**
 * Setup subsolver object using supplied cryptset
 *
 * @param wxc Cryptset to use for setting
 */
    void setup(const WxCryptSet& wxc, DM dm);

/**
 * Step the solver by given time step. If this step failed and if the
 * parent solver is running in variable time-stepping mode, it will be
 * called again with a smaller time-step. Otherwise, on failure the
 * simulation will be aborted.
 *
 * @param dt time step
 *
 * @return an instance of a crypt. This should contain the following
 * information about this step.
 * 'Status' : [required] bool. True if success, false if failure
 * 'Dt' : [optional] double. A suggestion for time-step on next step.
 * 'Message' : [optional] string. A message indicating why step failed (if it did)
 */
    virtual WxStepperStatus<REAL> step(REAL t, REAL dt, Vec in, Vec out);

/**
 * Return limited value.
 *
 */
    void applyToVector(wxNodalDGgeometry2D<REAL> *geom, Vec qk, Vec q_limited);

  protected:

/**
 * Apply Limiter to given vector
 *
 * @param dm [in] data management object
 * @param geom [in] DG geometry object
 * @param qk [in] conserved variables to which limiter is to be applied
 * @param q_limited [out] comverved variables after the limiter has been applied
 */
    virtual void applyLimiter(wxNodalDGgeometry2D<REAL> *geom, Vec qk, Vec q_limited) = 0;

  private:

};

#endif // WXNODALDGLIMITER_H
