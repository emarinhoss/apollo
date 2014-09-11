#ifndef APSOLVER_H
#define APSOLVER_H

// WarpX solver includes
#include <wxsolverbase.h>
#include <apsubsolver.h>
#include <apsubsolverstep.h>

// PETSc includes
#include <petscdmplex.h>

// std includes
#include <map>
#include <string>

typedef struct {
  Vec dg_vars;
  Vec cg_vars;
  char filename[PETSC_MAX_PATH_LEN];
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
    virtual WxStepperStatus<REAL> step(REAL dt);


    DM getdatamanagment(){
        return _dm;
    }

    UserContext getusercontext(){
        return _usr;
    }

    PetscViewer getViewer(){
        return _viewer;
    }

  private:

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
 *  output the data in vtk format
 *
 */
   void OutputVTK(DM dm, unsigned frame);

/**
 * Advance solution from tbeg to tend
 *
 * @param tbeg Starting time for advance
 * @param tend End time for advance
 * @param initial time-step to use. On return the last good time step
 */
    void advance(REAL tbeg, REAL tend, REAL& dt);

/**
 * Read a mesh and create data management object
 */
    void createMesh(MPI_Comm comm, DM *dm);

/**
 * check the grid
 */
    void SetupLocalSpace(DM dm, UserContext usr);

/**
 * Write comboSolver data to node
 *
 * @param io Pointer to file node
 * @param frame Frame number to write
 * @param tcurr Time at which data is written
 * @param telapsed Time for this advance
 */
//    void writeData(WxIoBase *io, unsigned frame, REAL tcurr, REAL telapsed);

/** Types of string -> WxSubSolver */
    typedef std::map<std::string, ApSubSolver<REAL>* > SubSolverMap_t;
    typedef std::pair<std::string, ApSubSolver<REAL>* > SubSolverPair_t;

/** Types of string -> WxSubSolverStep */
    typedef std::map<std::string, ApSubSolverStep<REAL> > SubSolverStepMap_t;
    typedef std::pair<std::string, ApSubSolverStep<REAL> > SubSolverStepPair_t;


/** Start and end time of simulation */
    REAL _tstart, _tend;
/** First frame number */
    unsigned _startFrame;
/** No of output files to write */
    unsigned _nout;
/** Initial time-step to use */
    REAL _dt;
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
/** Petsc Data Management object */
    DM _dm;
/** Petsc error code */
    PetscErrorCode  _ierr;
/** Grid filename to be read and name of output file*/
    char *_filename, *_outfname;
/** handle to visualize data */
    PetscViewer _viewer;
/** user-defined context */
    UserContext _usr;
/** domain partinioner */
    char *_partitioner;

};


#endif // APSOLVER_H
