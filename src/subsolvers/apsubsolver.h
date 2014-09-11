#ifndef __apsubsolver__h__
#define __apsubsolver__h__

// WarpX lib includes
#include <wxcrypt.h>
#include <wxmsgbase.h>
#include <wxobject.h>
#include <wxstepper.h>

// std includes
#include <string>
#include <typeinfo>
#include <vector>

// forward declare solver
template <typename REAL> class ApSolver;

/**
 * Base class for sub-solvers in WarpX system. A subsolver is an
 * abject which encodes algorithms to transform input variables into
 * output variables. It is like a subroutine in a programming
 * langauge.
 */
template <typename REAL>
class ApSubSolver : public WxStepper<REAL>
{
  public:
/**
 * Create new subsolver with given name
 *
 * @param name Name of subsolver
 */
    ApSubSolver(const std::string& name);

/** Destroy subsolver */
    virtual ~ApSubSolver();

/** Return type of subsolver */
    virtual const std::type_info& type() const;

/**
 * Setup subsolver object using supplied cryptset
 *
 * @param wxc Cryptset to use for setting
 */
    virtual void setup(const WxCryptSet& wxc);

/**
 * Set parent solver object
 * 
 * @param parent pointer to parent
 */
    void setParent(ApSolver<REAL> *parent);

/**
 * Get parent solver object
 * 
 * @return pointer to parent
 */
    ApSolver<REAL>* getParent() const;


  private:
    ApSolver<REAL> *_parent; // parent solver

};

#endif //  __apsubsolver__h__
