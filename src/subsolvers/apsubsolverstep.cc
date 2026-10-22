// WarpX subsolver includes
#include "apsubsolverstep.h"
#include "apsubsolver.h"
#include "apsolver.h"

// WarpX lib includes
#include <wxany.h>
#include <wxlogger.h>
#include <wxlogstream.h>

// std includes
#include <sstream>
#include <vector>
#include <string>
#include <iostream>

template <typename REAL>
void
ApSubSolverStep<REAL>::setup(const WxCryptSet& wxc, DM dm)
{
    WxLogger *l = WxLogger::get("apollo-root.console");
    WxLogStream infoStrm = l->getInfoStream();
    infoStrm << "Setting up subsolver-step '" << wxc.name() << "'" << std::endl;

    // time fraction
    if (wxc.has("DtFrac"))
        dtFrac = wxc.template get<REAL>("DtFrac");
    else
        dtFrac = 1.0;

    std::vector<WxAny>::const_iterator itr;
    // list of subsolvers to apply
    std::vector<WxAny> ssa = wxc.template get<std::vector<WxAny> >("SubSolvers");

    for (itr=ssa.begin(); itr!=ssa.end(); ++itr)
        subSolvers.push_back( wx_any_cast<std::string>(*itr) );

}

// instantiations
//template class ApSubSolverStep<float>;
template class ApSubSolverStep<double>;
