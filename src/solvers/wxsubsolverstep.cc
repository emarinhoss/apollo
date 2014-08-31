// WarpX subsolver includes
#include "wxsubsolverstep.h"
#include "wxsubsolver.h"
#include "apsolver.h"

// WarpX lib includes
#include <wxany.h>

// std includes
#include <sstream>
#include <vector>
#include <string>
#include <iostream>

template <typename REAL>
void
WxSubSolverStep<REAL>::setup(const WxCryptSet& wxc)
{
    //debStrm << "Setting up subsolver-step '" << wxc.name() << "'" << std::endl;

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
template class WxSubSolverStep<float>;
template class WxSubSolverStep<double>;
