// includes
#include <wxcreator.h>

// hyper includes
#include <apdomaindecompcheck.h>

//// hyperSubSolver
WxCreator<ApDomainDecompCheck<float>, ApSubSolver<float> > __domaindecomp_f("checkDomainDecomp");
WxCreator<ApDomainDecompCheck<double>, ApSubSolver<double> > __domaindecomp_d("checkDomainDecomp");
