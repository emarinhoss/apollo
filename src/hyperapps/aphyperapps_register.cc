// WarpX includes
#include <wxcreator.h>

// Advection Problems includes
#include "advection/wxadvectioneqn.h"
WxCreator< WxAdvectionEqn<float>, WxHyperbolicEqn<float> > __advectionEqn_f("advectionEqn");
WxCreator< WxAdvectionEqn<double>, WxHyperbolicEqn<double> > __advectionEqn_d("advectionEqn");

#include "advection/wxadvectionunstuctured.h"
WxCreator< WxAdvectionUnstructuredEqn<float>, WxHyperbolicEqn<float> > __advectionUnsEqn_f("advectionUnstEqn");
WxCreator< WxAdvectionUnstructuredEqn<double>, WxHyperbolicEqn<double> > __advectionUnsEqn_d("advectionUnstEqn");


// Euler Equation
//#include "euler/wmeulereqn.h"
//WxCreator< WmEulerEqn<float>, WmHyperbolicEqn<float> > __eulerEqn_f("eulerEqn_f");
//WxCreator< WmEulerEqn<double>, WmHyperbolicEqn<double> > __eulerEqn_d("eulerEqn_d");

//WxCreator< WmEulerEqn<float>::WmEulerBC_Inflow,     WmEquationBC<float> > __eulerEqnInflow_f("eulerEqnBC_Inflow_f");
//WxCreator< WmEulerEqn<float>::WmEulerBC_Outflow,    WmEquationBC<float> > __eulerEqnOutflow_f("eulerEqnBC_Outflow_f");
//WxCreator< WmEulerEqn<float>::WmEulerBC_SolidWall,  WmEquationBC<float> > __eulerEqnSolidWall_f("eulerEqnBC_SolidWall_f");
//WxCreator< WmEulerEqn<double>::WmEulerBC_Inflow,    WmEquationBC<double> > __eulerEqnInflow_d("eulerEqnBC_Inflow_d");
//WxCreator< WmEulerEqn<double>::WmEulerBC_Outflow,   WmEquationBC<double> > __eulerEqnOutflow_d("eulerEqnBC_Outflow_d");
//WxCreator< WmEulerEqn<double>::WmEulerBC_SolidWall, WmEquationBC<double> > __eulerEqnSolidWall_d("eulerEqnBC_SolidWall_d");
