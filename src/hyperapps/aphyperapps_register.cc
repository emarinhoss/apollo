// WarpX includes
#include <wxcreator.h>

// Advection Problems includes
#include "advection/wxadvectioneqn.h"
//WxCreator< WxAdvectionEqn<float>, WxHyperbolicEqn<float> > __advectionEqn_f("advectionEqn");
WxCreator< WxAdvectionEqn<double>, WxHyperbolicEqn<double> > __advectionEqn_d("advectionEqn");

#include "advection/wxadvectionunstuctured.h"
//WxCreator< WxAdvectionUnstructuredEqn<float>, WxHyperbolicEqn<float> > __advectionUnsEqn_f("advectionUnstEqn");
WxCreator< WxAdvectionUnstructuredEqn<double>, WxHyperbolicEqn<double> > __advectionUnsEqn_d("advectionUnstEqn");

#include "advection/wxgridfixedbc.h"
//WxCreator< WxGridFixedBC<float>, ApSubSolver<float> > __fixedVBC_f("gridFixedBC");
WxCreator< WxGridFixedBC<double>, ApSubSolver<double> > __fixedVBC_d("gridFixedBC");

/**
 * Euler
 */

// ------------------ Equations --------------------------//
#include "euler/wxeulereqn.h"
//WxCreator< WxEulerEqn<float>, WxHyperbolicEqn<float> > __eulerEqn_f("eulerEqn");
WxCreator< WxEulerEqn<double>, WxHyperbolicEqn<double> > __eulerEqn_d("eulerEqn");

#include "euler/apeulerentropyeqn.h"
//WxCreator< WxEulerEntropyEqn<float>, WxHyperbolicEqn<float> > __eulerEntEqn_f("eulerEntropyEqn");
WxCreator< WxEulerEntropyEqn<double>, WxHyperbolicEqn<double> > __eulerEntEqn_d("eulerEntropyEqn");

// ------------------ Boundary Conditions --------------------------//
#include "euler/wxisentropicvortexbc.h"
//WxCreator< WxIsentropicVortexBC<float>, ApSubSolver<float> > __isentropicVBC_f("isentropicVortexBC");
WxCreator< WxIsentropicVortexBC<double>, ApSubSolver<double> > __isentropicVBC_d("isentropicVortexBC");

#include "euler/wxeulerdirichletbc.h"
//WxCreator< WxEulerDirichletBC<float>, ApSubSolver<float> > __edirichletBC_f("eulerDirichletBC");
WxCreator< WxEulerDirichletBC<double>, ApSubSolver<double> > __edirichletBC_d("eulerDirichletBC");

#include "euler/wxeulerinflowbc.h"
//WxCreator< WxEulerInflowBC<float>, ApSubSolver<float> > __eInflowBC_f("eulerInflowBC");
WxCreator< WxEulerInflowBC<double>, ApSubSolver<double> > __eInflowBC_d("eulerInflowBC");

#include "euler/wxeulerwallbc.h"
//WxCreator< WxEulerWallBC<float>, ApSubSolver<float> > __eWallBC_f("eulerWallBC");
WxCreator< WxEulerWallBC<double>, ApSubSolver<double> > __eWallBC_d("eulerWallBC");

#include "euler/wxeulerzerogradientbc.h"
//WxCreator< WxEulerZeroGradientBC<float>, ApSubSolver<float> > __eZGBC_f("eulerZeroGradientBC");
WxCreator< WxEulerZeroGradientBC<double>, ApSubSolver<double> > __eZGBC_d("eulerZeroGradientBC");

#include "euler/apeulerinflow_entropybc.h"
//WxCreator< WxEulerInflowEntropyBC<float>, ApSubSolver<float> > __eIEBC_f("eulerInflowEntropyBC");
WxCreator< WxEulerInflowEntropyBC<double>, ApSubSolver<double> > __eIEBC_d("eulerInflowEntropyBC");

// ------------------ Limiters --------------------------//
#include "euler/wxhesthavenwarburtoneulerlimiter.h"
//WxCreator< WxHestavenWarburtonEulerLimiter<float>, ApSubSolver<float> > __elimHW_f("eulerLimiterHW");
WxCreator< WxHestavenWarburtonEulerLimiter<double>, ApSubSolver<double> > __elimHW_d("eulerLimiterHW");

// ------------------ Source Terms --------------------------//
#include "euler/aprmfelectronvelocitysrc.h"
//WxCreator< ApRMFelectronVelocitySrc<float>, WxHyperbolicSrc<float> > __rmfevsrc_f("rmfElectronVelSrc");
WxCreator< ApRMFelectronVelocitySrc<double>, WxHyperbolicSrc<double> > __rmfevsrc_d("rmfElectronVelSrc");

/**
 * Maxwell
 */

// Equations
#include "maxwell/wxphmaxwelleqn.h"
//WxCreator< WxPHMaxwellEqn<float>, WxHyperbolicEqn<float> > __maxwellEqn_f("phMaxwellEqn");
WxCreator< WxPHMaxwellEqn<double>, WxHyperbolicEqn<double> > __maxwellEqn_d("phMaxwellEqn");

// ------------------ Boundary Conditions --------------------------//
#include "maxwell/wxmaxwelltranversemagneticbc.h"
//WxCreator< WxMaxwellTransverseMagnetic<float>, ApSubSolver<float> > __mTMBC_f("maxwellTransverseMagneticBC");
WxCreator< WxMaxwellTransverseMagnetic<double>, ApSubSolver<double> > __mTMBC_d("maxwellTransverseMagneticBC");

#include "maxwell/wxphmaxwellconductingwallbc.h"
//WxCreator< WxPHMaxwellConductingWallBC<float>, ApSubSolver<float> > __mCWBC_f("phmConductingWallBC");
WxCreator< WxPHMaxwellConductingWallBC<double>, ApSubSolver<double> > __mCWBC_d("phmConductingWallBC");

#include "maxwell/wxphmaxwellopenbc.h"
//WxCreator< WxPHMaxwellOpenBC<float>, ApSubSolver<float> > __phOBC_f("phmOpenBC");
WxCreator< WxPHMaxwellOpenBC<double>, ApSubSolver<double> > __phOBC_d("phmOpenBC");

#include "maxwell/wxmaxwellrmfbc.h"
//WxCreator< WxMaxwellRMFBC<float>, ApSubSolver<float> > __phmRMFBC_f("maxwellRMFBC");
WxCreator< WxMaxwellRMFBC<double>, ApSubSolver<double> > __phmRMFBC_d("maxwellRMFBC");

#include "maxwell/apmaxwellrmfantennabc.h"
//WxCreator< WxMaxwellRMFAntennaBC<float>, ApSubSolver<float> > __phmRMFantBC_f("maxwellRMFAntennaBC");
WxCreator< WxMaxwellRMFAntennaBC<double>, ApSubSolver<double> > __phmRMFantBC_d("maxwellRMFAntennaBC");

#include "maxwell/apphmaxwelldielectricbc.h"
//WxCreator< WxPHMaxwellDielectricBC<float>, ApSubSolver<float> > __phmDieBC_f("phmDielectricBC");
WxCreator< WxPHMaxwellDielectricBC<double>, ApSubSolver<double> > __phmDieBC_d("phmDielectricBC");

#include "maxwell/apemfielsrmfwithharmonics.h"
//WxCreator< ApEMFieldsRMFWithHarmonics<float>, ApSubSolver<float> > __phmrmfWHBC_f("maxwellRMFHarmonicsBC");
WxCreator< ApEMFieldsRMFWithHarmonics<double>, ApSubSolver<double> > __phmrmfWHBC_d("maxwellRMFHarmonicsBC");

#include "maxwell/apsimplifiedrmfbc.h"
//WxCreator< APsimplifiedRMFbc<float>, ApSubSolver<float> > __phmsrmfBC_f("simplifiedRMFBC");
WxCreator< APsimplifiedRMFbc<double>, ApSubSolver<double> > __phmsrmfBC_d("simplifiedRMFBC");

// ------------------ Source Terms --------------------------//
#include "maxwell/wxpulsetrainsrc.h"
//WxCreator< WxPulseTrainSrc<float>, WxHyperbolicSrc<float> > __ptsrc_f("pulseTrain");
WxCreator< WxPulseTrainSrc<double>, WxHyperbolicSrc<double> > __ptsrc_d("pulseTrain");

#include "maxwell/wxtransversemagneticsrc.h"
//WxCreator< WxTransverseMagneticSrc<float>, WxHyperbolicSrc<float> > __tmsrc_f("tranverseMagneticSrc");
WxCreator< WxTransverseMagneticSrc<double>, WxHyperbolicSrc<double> > __tmsrc_d("tranverseMagneticSrc");

#include "maxwell/wxpointcurrent.h"
//WxCreator< WxpointCurrentSrc<float>, WxHyperbolicSrc<float> > __pcsrc_f("pointCurrentSrc");
WxCreator< WxpointCurrentSrc<double>, WxHyperbolicSrc<double> > __pcsrc_d("pointCurrentSrc");

#include "maxwell/aprmfsrc.h"
//WxCreator< ApRMFSrc<float>, WxHyperbolicSrc<float> > __aprmfsrc_f("maxwellRMFSrc");
WxCreator< ApRMFSrc<double>, WxHyperbolicSrc<double> > __aprmfsrc_d("maxwellRMFSrc");

#include "maxwell/aprmfantennasrc.h"
//WxCreator< ApRMFAntennaSrc<float>, WxHyperbolicSrc<float> > __aprmfant_f("maxwellRMFAntenna");
WxCreator< ApRMFAntennaSrc<double>, WxHyperbolicSrc<double> > __aprmfant_d("maxwellRMFAntenna");

/**
 * Multifluid sources, initializations, boundary conditions, and more.
 */

// ------------------ Source Terms --------------------------//
#include "multifluid/wxchargesrc.h"
//WxCreator< WxChargeSrc<float>, WxHyperbolicSrc<float> > __chargesrc_f("charges");
WxCreator< WxChargeSrc<double>, WxHyperbolicSrc<double> > __chargesrc_d("charges");

#include "multifluid/wxcurrentsrc.h"
//WxCreator< WxCurrentSrc<float>, WxHyperbolicSrc<float> > __currsrc_f("currents");
WxCreator< WxCurrentSrc<double>, WxHyperbolicSrc<double> > __currsrc_d("currents");

#include "multifluid/wxlorentzforcesrc.h"
//WxCreator< WxLorentzForceSrc<float>, WxHyperbolicSrc<float> > __lorForsrc_f("lorentzForces");
WxCreator< WxLorentzForceSrc<double>, WxHyperbolicSrc<double> > __lorForsrc_d("lorentzForces");

#include "multifluid/wxbraginskiifrictionsrc.h"
//WxCreator< WxBragFrictionSrc<float>, WxHyperbolicSrc<float> > __BFricSrc_f("bragFriction");
WxCreator< WxBragFrictionSrc<double>, WxHyperbolicSrc<double> > __BFricSrc_d("bragFriction");

#include "multifluid/apmagneticfluxcalc.h"
//WxCreator< ApMagneticFluxCalc<float>, ApAreaIntegral<float> > __MagFluxC_f("MagneticFlux");
WxCreator< ApMagneticFluxCalc<double>, ApAreaIntegral<double> > __MagFluxC_d("MagneticFlux");

#include "multifluid/ap2dmomentumxfer.h"
//WxCreator< Ap2DMomentumXfer<float>, WxHyperbolicSrc<float> > __2DMXfer_f("MomentumXfer2D");
WxCreator< Ap2DMomentumXfer<double>, WxHyperbolicSrc<double> > __2DMXfer_d("MomentumXfer2D");

#include "multifluid/apconstantresistivity.h"
//WxCreator< ApConstantResistivity<float>, WxHyperbolicSrc<float> > __ConsRes_f("constantResistivity");
WxCreator< ApConstantResistivity<double>, WxHyperbolicSrc<double> > __ConsRes_d("constantResistivity");

#include "multifluid/apanisotropicresistivitysrc.h"
//WxCreator< ApAnisotropicResistivitySrc<float>, WxHyperbolicSrc<float> > __anisoRes_f("anisotropicResistivity");
WxCreator< ApAnisotropicResistivitySrc<double>, WxHyperbolicSrc<double> > __anisoRes_d("anisotropicResistivity");

#include "multifluid/apcollisionalradiativemodelingsrc.h"
//WxCreator< ApCollisionalRadiativeModelingSrc<float>, WxHyperbolicSrc<float> > __colRadMod_f("CRmodel");
WxCreator< ApCollisionalRadiativeModelingSrc<double>, WxHyperbolicSrc<double> > __colRadMod_d("CRmodel");

// ------------------ Boundary Conditions --------------------------//
#include "multifluid/wxtwofluidrmfbc.h"
//WxCreator< WxTwoFluidRMFBC<float>, ApSubSolver<float> > __tfrmfBC_f("twoFluidRMFBC");
WxCreator< WxTwoFluidRMFBC<double>, ApSubSolver<double> > __tfrmfBC_d("twoFluidRMFBC");

#include "multifluid/wxtwofluidconductingwallbc.h"
//WxCreator< WxTwoFluidConductingWallBC<float>, ApSubSolver<float> > __tfcwBC_f("twoFluidConductingWallBC");
WxCreator< WxTwoFluidConductingWallBC<double>, ApSubSolver<double> > __tfcwBC_d("twoFluidConductingWallBC");

#include "multifluid/aptwofluiddielectricbc.h"
//WxCreator< WxTwoFluidDielectricBC<float>, ApSubSolver<float> > __tfdBC_f("twoFluidDielectricBC");
WxCreator< WxTwoFluidDielectricBC<double>, ApSubSolver<double> > __tfdBC_d("twoFluidDielectricBC");

#include "multifluid/aptwofluidrmfantennabc.h"
//WxCreator< WxTwoFluidRMFAntennaBC<float>, ApSubSolver<float> > __tfrmfaBC_f("twoFluidAntennaBC");
WxCreator< WxTwoFluidRMFAntennaBC<double>, ApSubSolver<double> > __tfrmfaBC_d("twoFluidAntennaBC");

#include "multifluid/aptwofluidrmfwithharmonics.h"
//WxCreator< WxTwoFluidRMFWithHarmonicsBC<float>, ApSubSolver<float> > __tfrmfwhBC_f("twoFluidRMFHarmonicsBC");
WxCreator< WxTwoFluidRMFWithHarmonicsBC<double>, ApSubSolver<double> > __tfrmfwhBC_d("twoFluidRMFHarmonicsBC");

#include "multifluid/aptwofluidsimplifiedrmfbc.h"
//WxCreator< APTwoFluidSimplifiedRMFBC<float>, ApSubSolver<float> > __tfsrmfBC_f("twoFluidSimplifiedRMFBC");
WxCreator< APTwoFluidSimplifiedRMFBC<double>, ApSubSolver<double> > __tfsrmfBC_d("twoFluidSimplifiedRMFBC");

#include "multifluid/aptwofluidomfbc.h"
//WxCreator< APTwoFluidOMFBC<float>, ApSubSolver<float> > __tfomfBC_f("twoFluidOMFBC");
WxCreator< APTwoFluidOMFBC<double>, ApSubSolver<double> > __tfomfBC_d("twoFluidOMFBC");

#include "multifluid/apthreefluidconstantbc.h"
//WxCreator< ApThreeFluidConstantBC<float>, ApSubSolver<float> > __tfcBC_f("threeFluidConstantBC");
WxCreator< ApThreeFluidConstantBC<double>, ApSubSolver<double> > __tfcBC_d("threeFluidConstantBC");

// ------------------ Limiters --------------------------//
#include "multifluid/aptualiabadilimiter.h"
//WxCreator< WxTuAliabadiLimiter<float>, ApSubSolver<float> > __taalim_f("tuAliabadiLimiter");
WxCreator< WxTuAliabadiLimiter<double>, ApSubSolver<double> > __taalim_d("tuAliabadiLimiter");
