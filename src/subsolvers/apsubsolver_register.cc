// includes
#include <wxcreator.h>

// hyper includes
#include "apdomaindecompcheck.h"
WxCreator<ApDomainDecompCheck<float>, ApSubSolver<float> > __domaindecomp_f("checkDomainDecomp");
WxCreator<ApDomainDecompCheck<double>, ApSubSolver<double> > __domaindecomp_d("checkDomainDecomp");

#include "apfvm2dscheme.h"
WxCreator<ApFVM2Dscheme<float>, ApSubSolver<float> > __FVM2D_f("fv2dscheme");
WxCreator<ApFVM2Dscheme<double>, ApSubSolver<double> > __FVM2D_d("fv2dscheme");
