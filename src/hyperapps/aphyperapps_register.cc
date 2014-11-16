// WarpX includes
#include <wxcreator.h>

// Advection Problems includes
#include "advection/wxadvectioneqn.h"
WxCreator< WxAdvectionEqn<float>, WxHyperbolicEqn<float> > __advectionEqn_f("advectionEqn");
WxCreator< WxAdvectionEqn<double>, WxHyperbolicEqn<double> > __advectionEqn_d("advectionEqn");

#include "advection/wxadvectionunstuctured.h"
WxCreator< WxAdvectionUnstructuredEqn<float>, WxHyperbolicEqn<float> > __advectionUnsEqn_f("advectionUnstEqn");
WxCreator< WxAdvectionUnstructuredEqn<double>, WxHyperbolicEqn<double> > __advectionUnsEqn_d("advectionUnstEqn");

#include "advection/wxgridfixedbc.h"
WxCreator< WxGridFixedBC<float>, WxGridBC<float> > __fixedVBC_f("gridFixedBC");
WxCreator< WxGridFixedBC<float>, WxGridBC<float> > __fixedVBC_d("gridFixedBC");

// Euler Equation
#include "euler/wxeulereqn.h"
WxCreator< WxEulerEqn<float>, WxHyperbolicEqn<float> > __eulerEqn_f("eulerEqn");
WxCreator< WxEulerEqn<double>, WxHyperbolicEqn<double> > __eulerEqn_d("eulerEqn");
