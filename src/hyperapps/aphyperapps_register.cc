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
WxCreator< WxGridFixedBC<float>, ApSubSolver<float> > __fixedVBC_f("gridFixedBC");
WxCreator< WxGridFixedBC<double>, ApSubSolver<double> > __fixedVBC_d("gridFixedBC");

/**
 * Euler
 */

// ------------------ Equations --------------------------//
#include "euler/wxeulereqn.h"
WxCreator< WxEulerEqn<float>, WxHyperbolicEqn<float> > __eulerEqn_f("eulerEqn");
WxCreator< WxEulerEqn<double>, WxHyperbolicEqn<double> > __eulerEqn_d("eulerEqn");

// ------------------ Boundary Conditions --------------------------//
#include "euler/wxisentropicvortexbc.h"
WxCreator< WxIsentropicVortexBC<float>, ApSubSolver<float> > __isentropicVBC_f("isentropicVortexBC");
WxCreator< WxIsentropicVortexBC<double>, ApSubSolver<double> > __isentropicVBC_d("isentropicVortexBC");

#include "euler/wxeulerdirichletbc.h"
WxCreator< WxEulerDirichletBC<float>, ApSubSolver<float> > __edirichletBC_f("eulerDirichletBC");
WxCreator< WxEulerDirichletBC<double>, ApSubSolver<double> > __edirichletBC_d("eulerDirichletBC");

#include "euler/wxeulerinflowbc.h"
WxCreator< WxEulerInflowBC<float>, ApSubSolver<float> > __eInflowBC_f("eulerInflowBC");
WxCreator< WxEulerInflowBC<double>, ApSubSolver<double> > __eInflowBC_d("eulerInflowBC");

#include "euler/wxeulerwallbc.h"
WxCreator< WxEulerWallBC<float>, ApSubSolver<float> > __eWallBC_f("eulerWallBC");
WxCreator< WxEulerWallBC<double>, ApSubSolver<double> > __eWallBC_d("eulerWallBC");

#include "euler/wxeulerzerogradientbc.h"
WxCreator< WxEulerZeroGradientBC<float>, ApSubSolver<float> > __eZGBC_f("eulerZeroGradientBC");
WxCreator< WxEulerZeroGradientBC<double>, ApSubSolver<double> > __eZGBC_d("eulerZeroGradientBC");

// ------------------ Limiters --------------------------//
#include "euler/wxhesthavenwarburtoneulerlimiter.h"
WxCreator< WxHestavenWarburtonEulerLimiter<float>, ApSubSolver<float> > __elimHW_f("eulerLimiterHW");
WxCreator< WxHestavenWarburtonEulerLimiter<double>, ApSubSolver<double> > __elimHW_d("eulerLimiterHW");

/**
 * Maxwell
 */

// Equations
#include "maxwell/wxphmaxwelleqn.h"
WxCreator< WxPHMaxwellEqn<float>, WxHyperbolicEqn<float> > __maxwellEqn_f("phMaxwellEqn");
WxCreator< WxPHMaxwellEqn<double>, WxHyperbolicEqn<double> > __maxwellEqn_d("phMaxwellEqn");

// ------------------ Boundary Conditions --------------------------//
#include "maxwell/wxmaxwelltranversemagneticbc.h"
WxCreator< WxMaxwellTransverseMagnetic<float>, ApSubSolver<float> > __mTMBC_f("maxwellTransverseMagneticBC");
WxCreator< WxMaxwellTransverseMagnetic<double>, ApSubSolver<double> > __mTMBC_d("maxwellTransverseMagneticBC");

#include "maxwell/wxphmaxwellconductingwallbc.h"
WxCreator< WxPHMaxwellConductingWallBC<float>, ApSubSolver<float> > __mCWBC_f("phmConductingWallBC");
WxCreator< WxPHMaxwellConductingWallBC<double>, ApSubSolver<double> > __mCWBC_d("phmConductingWallBC");

#include "maxwell/wxphmaxwellopenbc.h"
WxCreator< WxPHMaxwellOpenBC<float>, ApSubSolver<float> > __phOBC_f("phmOpenBC");
WxCreator< WxPHMaxwellOpenBC<double>, ApSubSolver<double> > __phOBC_d("phmOpenBC");

#include "maxwell/wxmaxwellrmfbc.h"
WxCreator< WxMaxwellRMFBC<float>, ApSubSolver<float> > __phmRMFBC_f("maxwellRMFBC");
WxCreator< WxMaxwellRMFBC<double>, ApSubSolver<double> > __phmRMFBC_d("maxwellRMFBC");

#include "maxwell/apmaxwellrmfantennabc.h"
WxCreator< WxMaxwellRMFAntennaBC<float>, ApSubSolver<float> > __phmRMFantBC_f("maxwellRMFAntennaBC");
WxCreator< WxMaxwellRMFAntennaBC<double>, ApSubSolver<double> > __phmRMFantBC_d("maxwellRMFAntennaBC");

// ------------------ Source Terms --------------------------//
#include "maxwell/wxpulsetrainsrc.h"
WxCreator< WxPulseTrainSrc<float>, WxHyperbolicSrc<float> > __ptsrc_f("pulseTrain");
WxCreator< WxPulseTrainSrc<double>, WxHyperbolicSrc<double> > __ptsrc_d("pulseTrain");

#include "maxwell/wxtransversemagneticsrc.h"
WxCreator< WxTransverseMagneticSrc<float>, WxHyperbolicSrc<float> > __tmsrc_f("tranverseMagneticSrc");
WxCreator< WxTransverseMagneticSrc<double>, WxHyperbolicSrc<double> > __tmsrc_d("tranverseMagneticSrc");

#include "maxwell/wxpointcurrent.h"
WxCreator< WxpointCurrentSrc<float>, WxHyperbolicSrc<float> > __pcsrc_f("pointCurrentSrc");
WxCreator< WxpointCurrentSrc<double>, WxHyperbolicSrc<double> > __pcsrc_d("pointCurrentSrc");


/**
 * Multifluid sources, initializations, boundary conditions, and more.
 */

// ------------------ Source Terms --------------------------//
#include "multifluid/wxchargesrc.h"
WxCreator< WxChargeSrc<float>, WxHyperbolicSrc<float> > __chargesrc_f("charges");
WxCreator< WxChargeSrc<double>, WxHyperbolicSrc<double> > __chargesrc_d("charges");

#include "multifluid/wxcurrentsrc.h"
WxCreator< WxCurrentSrc<float>, WxHyperbolicSrc<float> > __currsrc_f("currents");
WxCreator< WxCurrentSrc<double>, WxHyperbolicSrc<double> > __currsrc_d("currents");

#include "multifluid/wxlorentzforcesrc.h"
WxCreator< WxLorentzForceSrc<float>, WxHyperbolicSrc<float> > __lorForsrc_f("lorentzForces");
WxCreator< WxLorentzForceSrc<double>, WxHyperbolicSrc<double> > __lorForsrc_d("lorentzForces");

#include "multifluid/wxbraginskiifrictionsrc.h"
WxCreator< WxBragFrictionSrc<float>, WxHyperbolicSrc<float> > __BFricSrc_f("bragFriction");
WxCreator< WxBragFrictionSrc<double>, WxHyperbolicSrc<double> > __BFricSrc_d("bragFriction");

#include "multifluid/apmagneticfluxcalc.h"
WxCreator< ApMagneticFluxCalc<float>, ApAreaIntegral<float> > __MagFluxC_f("MagneticFlux");
WxCreator< ApMagneticFluxCalc<double>, ApAreaIntegral<double> > __MagFluxC_d("MagneticFlux");

// ------------------ Boundary Conditions --------------------------//
#include "multifluid/wxtwofluidrmfbc.h"
WxCreator< WxTwoFluidRMFBC<float>, ApSubSolver<float> > __tfrmfBC_f("twoFluidRMFBC");
WxCreator< WxTwoFluidRMFBC<double>, ApSubSolver<double> > __tfrmfBC_d("twoFluidRMFBC");

#include "multifluid/wxtwofluidconductingwallbc.h"
WxCreator< WxTwoFluidConductingWallBC<float>, ApSubSolver<float> > __tfcwBC_f("twoFluidConductingWallBC");
WxCreator< WxTwoFluidConductingWallBC<double>, ApSubSolver<double> > __tfcwBC_d("twoFluidConductingWallBC");
