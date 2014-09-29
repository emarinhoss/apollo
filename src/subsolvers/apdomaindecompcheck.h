#ifndef APDOMAINDECOMPCHECK_H
#define APDOMAINDECOMPCHECK_H

// WarpX lib includes
#include <wxfunction.h>

// WarpX subsolver includes
#include "apsubsolver.h"
#include <apsolver.h>

// std includes
#include <vector>


/**
 * Computes central difference of read arrays
 */
template <typename REAL>
class ApDomainDecompCheck : public ApSubSolver<REAL>
{
  public:
    ApDomainDecompCheck()
            : ApSubSolver<REAL>("checkDomainDecomp") {
    }

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
    void init();

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

  private:
    DM _dm;
    UserContext _usr;


};
#endif // APDOMAINDECOMPCHECK_H
