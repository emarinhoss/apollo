#ifndef __plsimulation__
#define __plsimulation__

// lib includes
#include <wxobject.h>
#include <apsolver.h>

// Mpi includes
# include <mpi.h>
# include <wxmpimsg.h>
#include <petsc.h>

/**
 * Top level class for Apollo simulations. This class does all the
 * setup of the simulation and runs it to completion.
 */
template <typename REAL>
class ApSimulation : public WxObject
{
  public:
/** Construct a new simulation object */
    ApSimulation(int argc, char **argv);

/** Destroy simulation */
    virtual ~ApSimulation();

/**
 * Return name of input file
 *
 * @return Name of input file
 */
    std::string getInpFileName() const;

/**
 * Return name of run
 *
 * @return Name of run
 */
    std::string getRunName() const;

/**
 * Return true if we are restarting simulation
 *
 * @return true if we are restarting simulation
 */
    bool isRestarting() const;

/**
 * Get file from which to restart the simulation.
 *
 * @param file name
 */
    std::string getRestartFile() const;

/**
 * Setup simulation using supplied crypset.
 *
 * @param wxc Cryptset using which the object is set up.
 */
    virtual void setup(const WxCryptSet& wxc);

/**
 * Initialize simulation.
 */
    virtual void init();

/**
 * Load simulation from file.
 *
 * @param io I/O object to use for reading
 * @param grpNode group node to read from
 */
//    virtual void load(WxIoBase& io, const WxIoNodeType& grpNode);

/** Run simulation */
    void simulate();

  private:
/** Top level solver in simulation */
    ApSolver<REAL> *solver;
/** Name of input file */
    std::string inpFileName;
/** Mesh file name */
    char meshFileName[PETSC_MAX_PATH_LEN];
/** Name of the run */
    std::string runName;
/** Set to true if output prefix is set */
    bool outPrefixSet;
/** Output prefix */
    std::string outPrefix;
/** Flag to indicate if simulation is a restart */
    bool restartSim;
/** File from which to restart */
    std::string restartFile;
/** Real number type to use */
    std::string realType;

/** Parse command line parameters */
    void parseCmdLine(int argc, char **argv);

/** Print usage message */
    void usage();

/** Remove characters following the last '.' in nm */
    std::string stripName(const std::string& nm);

    PetscInt _np; // number of processorszz
};

#endif // __plsimulation__
