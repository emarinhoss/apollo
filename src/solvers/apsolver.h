#ifndef APSOLVER_H
#define APSOLVER_H

// WarpX solver includes
#include "wxsolverbase.h"
#include <apsubsolver.h>
#include <apsubsolverstep.h>

// PETSc includes
#include <petscdmplex.h>
#include <petscts.h>
#include <wxpetsctimestepping.h>
#include <petscviewerhdf5.h>

// std includes
#include <map>
#include <string>

typedef struct {
  Vec dg_vars, cg_vars, solution;
  char filename[PETSC_MAX_PATH_LEN];
  unsigned nout;
  std::string runName;
} UserContext;

template <typename REAL>
class ApSolver : public WxSolverBase<REAL>
{
  public:
/**
 * Create solver with given name
 *
 * @param name Name of solver
 */
    ApSolver(const std::string& name);

/**
 * Destroy solver
 */
    virtual ~ApSolver();

/**
 * Set frame from which to start the simulation
 *
 * @param frame Frame from which to start
 */
    void setStartFrame(unsigned frame);

/**
 * Setup the subsolver using cryptset
 *
 * @param wxc Cryptset using which the object is set up.
 */
    virtual void setup(const WxCryptSet& wxc);

/**
 * Return pointer to subsolver
 *
 * @param name name of subsolver to return
 */
    ApSubSolver<REAL>* getSubSolver(const std::string& name);

/**
 * Get initial time-step to use
 *
 * @return initial time-step to use
 */
    REAL getInitDt() const;

/**
 * Initialize the object. This is called after the setup() method is
 * called. This method is not called if the object is being restored
 * from an output file. In this case the load() method is called
 * instead.
 */
    virtual void init();

/**
 * Run solver.
 */
    virtual void solve();

/**
 * Advance the stepper by given time step. If this step failed and if
 * the parent solver is running in variable time-stepping mode, it
 * will be called again with a smaller time-step. A time-step
 * suggestion must be provided. The status of the step should be
 * encoded in the return object. See WxStepperStatus for details on
 * what fields can be set.
 *
 * @param dt time step to advance by
 * @return Status of stepper
 */
    virtual WxStepperStatus<REAL> step(REAL t, REAL dt, Vec in, Vec out);


    DM getdatamanagment(){
        return _dm;
    }

    UserContext getusercontext(){
        return _usr;
    }

//    PetscViewer getViewer(){
//        return _viewer;
//    }

    std::string getFilename_OutputVTK(){
        return _outfname;
    }

    /**
     *  output the data in vtk format
     *
     */
       void OutputVTK(DM dm, char *filename, PetscViewer *viewer);

/**
 * Writes a VTU file with the data
 *
 * @param [in] X - the vector to write to the .vtu file
 */
       void writeData(Vec X);

/**
 * Write a restart checkpoint for the current state.
 *
 * Called once per output frame, alongside the .vtu. The file is ROLLING - the
 * same name every time - so the whole mechanism costs one solution vector on
 * disk (5.2 MB on the phase3 mesh) rather than one per frame.
 *
 * @param X solution vector to checkpoint
 * @param time simulation time X is at
 */
       void writeCheckpoint(Vec X, REAL time);

/**
 * Restore state from a checkpoint written by writeCheckpoint().
 *
 * Call AFTER init(), not instead of it: init() creates the solution vector,
 * lays down the initial condition and sizes _dt from it, and the timestep must
 * come from the same place it came from on the original run or the resumed run
 * integrates differently from the one it claims to continue.
 *
 * Refuses, rather than resuming into a mismatch, when the checkpoint was
 * written by a different problem.
 *
 * @param path checkpoint file to read
 */
       void loadCheckpoint(const std::string& path);

/**
 * function that is to be used at every timestep
 * to display the iteration's progress.
 *
 * @param ts	- the TS context
 * @param steps- iteration number (after the final time step the
 *          monitor routine is called with a step of -1,
 *          this is at the final time which may have been interpolated to)
 * @param time	- current time
 * @param u	- current solution iterate
 * @param ctx	- [optional] monitoring context
 */
    PetscErrorCode MonitorVTK(TS ts, PetscInt stepnum, PetscReal time, Vec X, void *ctx);


/**
 * Sets the routine for evaluating the function, where U_t = F(t,u).
 *@param t	- current timestep
 *@param u	- input vector
 *@param F	- function vector
 *@param ctx	- [optional] user-defined function context
 */
    PetscErrorCode ComputeRHSforTS(TS ts, PetscReal t, Vec u, Vec global_out, void *ctx);

  private:

/**
 * Create a section that has all the data structure
 */
    void SetupLocalSpace(DM *dm);

/**
 * Type check each subsolver to ensure they have the correct variables
 * specified
 */
//    void typeCheck();

/**
 * Run startOnly steps
 */
    void startOnly();

/**
 * Run writeOnly steps
 *
 * @param io Pointer to I/O object
 */
//    void writeOnly(WxIoBase *io, WxIoNodeType ioNode);

/**
 * Run endOnly steps
 */
//    void endOnly();

/**
 * Read a mesh and create data management object
 */
    void createMesh(MPI_Comm comm, DM *dm);

/**
 * Set and get the output file name
 */
    void setFilename_OutputVTK(std::string fname){
        _outfname = fname;
    }

/** Types of string -> WxSubSolver */
    typedef std::map<std::string, ApSubSolver<REAL>* > SubSolverMap_t;
    typedef std::pair<std::string, ApSubSolver<REAL>* > SubSolverPair_t;

/** Types of string -> WxSubSolverStep */
    typedef std::map<std::string, ApSubSolverStep<REAL> > SubSolverStepMap_t;
    typedef std::pair<std::string, ApSubSolverStep<REAL> > SubSolverStepPair_t;


/** Start and end time of simulation */
    REAL _tstart, _tend, _tend_temp;
/** First frame number */
    unsigned _startFrame;
/** Current frame number */
    int _frameNum;
/** No of output intervals REMAINING to run. A restart shortens this. */
    unsigned _nout;
/** No of output files the DECK asked for. Set once, never changed.
 *
 * Kept apart from _nout because loadCheckpoint() shortens _nout to the
 * intervals that are left, and the checkpoint has to record what the deck
 * asked for. Writing the running value made restart non-idempotent: resuming
 * at frame 2 of 6 recorded "nout 4", and resuming from THAT checkpoint then
 * compared 4 against the deck's 6 and refused. One resume worked and a second
 * did not, which is the wrong half to get working. */
    unsigned _noutDeck;
/** Initial time-step to use */
    REAL _dt, _dt_temp;
/** Problem dimensions */
    PetscInt _dim;
/** Flag whether to use fuzzy stepper */
    int _useFixedDt;
/** Subsolvers used */
    SubSolverMap_t _subSolvers;
/** Steps to apply at start of simulation */
    std::vector<ApSubSolverStep<REAL> > _startOnly;
/** Steps to apply at end of simulation */
    std::vector<ApSubSolverStep<REAL> > _endOnly;
/** Steps to apply before writing */
    std::vector<ApSubSolverStep<REAL> > _writeOnly;
/** Sequence of steps for each time step */
    std::vector<ApSubSolverStep<REAL> > _perStep;

/** Petsc objects */
    DM _dm; // data management
    Vec solution;
    char *_filename; // gridfile name
    std::string _outfname; // output file name
//    PetscViewer _viewer; // viewer for output data
    UserContext _usr; // user-defined context
    int _fieldsNum; // Number of fields
    std::vector<std::string> _fieldsName;
    std::vector<int> _fieldsComponents;
    std::vector<int> _fieldsNumber;
    WxPetscTimeSteppingSolver<REAL, ApSolver> *tssolver;

};


#endif // APSOLVER_H
