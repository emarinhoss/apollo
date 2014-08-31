#ifndef __wxsolverbase__
#define __wxsolverbase__

// WarpX solver includes
#include <wxstepper.h>

// std includes
#include <string>

/**
 * A base class for solvers in WarpX.
 */
template <typename REAL>
class WxSolverBase : public WxStepper<REAL>
{
  public:
/**
 * Create new solver with given name
 *
 * @param name Name of solver
 */   
    WxSolverBase(const std::string& name);

/** 
 * Destroy the object
 */
    virtual ~WxSolverBase();

/**
 * Set the run name
 */
    void setRunName(const std::string& runName);

/**
 * Return the run name
 */
    std::string runName() const;

/**
 * Set the solver name
 *
 * @param name Solver name
 */
    void setSolverName(const std::string& name);

/**
 * Get the solver name
 *
 * @return Solver name
 */
    std::string getSolverName() const;

/**
 * Run the solver. This must be provided by child classes.
 */
    virtual void solve() = 0;

  private:
/** Name of the run. Used to create output files and logs */
    std::string _runName;
/** top solver name name */
    std::string _solverName;
};

#endif // __wxsolverbase__
