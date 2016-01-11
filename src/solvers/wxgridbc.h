#ifndef __wxgridbc__
#define __wxgridbc__

// WarpX subsolver includes
#include <apsubsolver.h>
#include "apsolver.h"
#include <wxpdggeometry.h>

/**
 * Base class which applies boundary conditions to fields living on a
 * grid box.
 */
template <typename REAL>
class WxGridBC : public ApSubSolver<REAL>
{
  public:

/**
 * Construct a new grid-bc object
 */
    WxGridBC(const std::string& name)
      : ApSubSolver<REAL>(name) {
    }

    virtual ~WxGridBC() {
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
 * Apply bounday condition to suppllied array. This is provided so one
 * can call BCs from inside other subsolvers
 *
 */    
    void applyToArray(REAL *xc, REAL *nx, REAL *q, REAL *qaux, REAL *AreaInts, REAL *qBC);

  protected:

/**
 * Apply BC to the lower edge along direction 'dir'
 *
 * @param dir direction in which to apply BC
 * @param arr array to which apply BC
 */
    virtual void applyBC(REAL *xc, REAL *nx, REAL *q, REAL *qaux, REAL *AreaInts, REAL *qBC) = 0;

  private:

};

#endif //  __wxgridbc__
