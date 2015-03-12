// includes
#include <wxcreator.h>

// hyper includes
#include "apdomaindecompcheck.h"
WxCreator<ApDomainDecompCheck<float>, ApSubSolver<float> > __domaindecomp_f("checkDomainDecomp");
WxCreator<ApDomainDecompCheck<double>, ApSubSolver<double> > __domaindecomp_d("checkDomainDecomp");

//#include "apfvm2dscheme.h"
//WxCreator<ApFVM2Dscheme<float>, ApSubSolver<float> > __FVM2D_f("fv2dscheme");
//WxCreator<ApFVM2Dscheme<double>, ApSubSolver<double> > __FVM2D_d("fv2dscheme");

#include "wxpgd2dscheme.h"
WxCreator<WxpDG2Dscheme<float>, ApSubSolver<float> > __DGM2D_f("dg2dscheme");
WxCreator<WxpDG2Dscheme<double>, ApSubSolver<double> > __DGM2D_d("dg2dscheme");

#include "wxnodaldg2dmethod.h"
WxCreator<WxNodalDG2dMethod<float>, ApSubSolver<float> > __FVM2D_f("nodalDG2d");
WxCreator<WxNodalDG2dMethod<double>, ApSubSolver<double> > __FVM2D_d("nodalDG2d");
