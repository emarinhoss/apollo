#ifndef WXPETSCTIMESTEPPING_H
#define WXPETSCTIMESTEPPING_H

//  petsc includes
#include <petscts.h>

// warpx includes
#include <wxlogger.h>
#include <wxlogstream.h>
#include <wxtimer.h>

/*
     User-defined application context
  */
template <typename REAL, typename CLS>
class WxPetscTimeSteppingSolver {
  public:
/**
 * Create a new SNES solver object.
 *
 * @param n Number of unknowns in equation system
 * @param cls Pointer to class providing functionality
 */

    WxPetscTimeSteppingSolver(DM dm, CLS *cls, MPI_Comm comm, REAL tstart, REAL tend, REAL dt)
    : dm(dm), cls(cls) {

      WxLogStream debStrm = WxLogger::get("warpx-root.console")->getDebugStream();

      // create time stepping scheme
      TSCreate(PETSC_COMM_WORLD, &solver);
      TSSetType(solver, TSSSP);

      TSMonitorSet(solver,WxPetscTimeSteppingSolver::MonitorVTK, (void*) cls,NULL);

      // relate the timestepping scheme
      // with the data managenent object
      TSSetDM(solver, dm);
      TSSetRHSFunction(solver,NULL,WxPetscTimeSteppingSolver::ComputeRHSforTS,(void*) cls);

      TSSetDuration(solver,1000,tend);
      TSSetInitialTimeStep(solver,tstart,dt);

    }

/**
 * Clean up memory when object is deleted
 */
    virtual ~WxPetscTimeSteppingSolver() {
      // allows output of jacobian in Matlab format
//       PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_MATLAB);
//       MatView(jacobian, PETSC_VIEWER_STDOUT_WORLD);

      TSDestroy(&solver);
      DMDestroy(&dm);
    }

/**
 * Find the roots of the the system of equations F(x) = b.
 *
 * @param qin Input guess for roots
 * @param qout Output solution
 * @return flag indicating if solution converged
 */
    bool solve(Vec X) {

        WxLogger *log = WxLogger::get("apollo-root.console");
        WxLogStream debStrm = log->getDebugStream();
        WxLogStream infStrm = log->getInfoStream();

        //REAL tcurr = this->getCurrentTime(); // starting time
        PetscScalar ftime;
        PetscInt nsteps;
        TSConvergedReason reason;

        TSSetSolution(solver,X);
        TSSetFromOptions(solver);
        TSSolve(solver,X);
        TSGetSolveTime(solver,&ftime);
        TSGetTimeStepNumber(solver,&nsteps);
        TSGetConvergedReason(solver,&reason);

        infStrm << TSConvergedReasons[reason] << " at time "
                << ftime << " after "
                << nsteps << " steps"
                << std::endl;

        return true;

    }

  private:

/**
 * Static method for evaluating function whose root is to be
 * found. Real work is done in evaluateFunction() method.
 *
 * @param snes SNES context
 * @param x Input vector.
 * @param f Output function f(x)
 * @param ctx Context for this method.
 * @return error code.
 */
    static PetscErrorCode MonitorVTK(TS ts, PetscInt stepnum, PetscReal time, Vec X, void *ctx) {
      return static_cast<CLS*>(ctx)->MonitorVTK(ts, stepnum, time, X, &ctx);
    }

/**
 * Static method for evaluating function whose root is to be
 * found. Real work is done in evaluateFunction() method.
 *
 * @param snes SNES context
 * @param x Input vector.
 * @param f Output function f(x)
 * @param ctx Context for this method.
 * @return error code.
 */
    static PetscErrorCode ComputeRHSforTS(TS ts,PetscReal t,Vec u,Vec F,void *ctx) {
      return static_cast<CLS*>(ctx)->ComputeRHSforTS(ts, t, u, F, &ctx);
    }

/** Pointer to service class */
    CLS *cls;
/** Data Management */
    DM dm;
/** Solver context */
    TS solver;

};

#endif // WXPETSCTIMESTEPPING_H
