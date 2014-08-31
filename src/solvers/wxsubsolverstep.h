#ifndef __wxsubsolverstep__h__
#define __wxsubsolverstep__h__

// WarpX includes
#include <wxobject.h>
#include <wxcryptset.h>

// std includes
#include <vector>
#include <string>

// forward declare WxComboSolver
template <typename REAL> class ApSolver;

template <typename REAL>
struct WxSubSolverStep : public WxObject
{

/**
 * Setup the subsolverstep using cryptset
 *
 * @param wxc Cryptset using which the object is set up.
 */
    virtual void setup(const WxCryptSet& wxc);

/** Fraction of time-step to apply */
    REAL dtFrac;
/** Subsolvers to apply */
    std::vector<std::string> subSolvers;
/** Variables to synchronize at end of step */
    std::vector<std::string> syncVars;
};

#endif //  __wxsubsolverstep__h__
