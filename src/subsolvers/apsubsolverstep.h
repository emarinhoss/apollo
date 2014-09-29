#ifndef __apsubsolverstep__h__
#define __apsubsolverstep__h__

// WarpX includes
#include <wxobject.h>
#include <wxcryptset.h>

// std includes
#include <vector>
#include <string>

// forward declare WxComboSolver
template <typename REAL> class ApSolver;

template <typename REAL>
struct ApSubSolverStep : public WxObject
{

/**
 * Setup the subsolverstep using cryptset
 *
 * @param wxc Cryptset using which the object is set up.
 */
    virtual void setup(const WxCryptSet& wxc, DM dm);

/** Fraction of time-step to apply */
    REAL dtFrac;
/** Subsolvers to apply */
    std::vector<std::string> subSolvers;
};

#endif //  __apsubsolverstep__h__
