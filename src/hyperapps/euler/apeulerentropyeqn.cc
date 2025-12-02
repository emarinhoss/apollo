#include "apeulerentropyeqn.h"

// WarpX lib includes
#include <wxlogger.h>

// WarpX hyperbolic solver includes
#include "wxeulereqn.h"
#include <wxlogger.h>
#include <wxlogstream.h>
// std includes
#include <cmath>

// flags to indicate which numerical flux to use
static const unsigned LF   = 0;
static const unsigned HLL  = 1;
static const unsigned ROE  = 2;
static const unsigned HLLC = 3;
static const unsigned WAVE = 4;

template<typename REAL>
void
WxEulerEntropyEqn<REAL>::
setup(const WxCryptSet& wxc)
{
    WxLogStream infoStrm = WxLogger::get("apollo-root.console")->getInfoStream();

  // set gas gamma
  _gas_gamma = wxc.template get<REAL>("gas_gamma");
  _efix = true;
  if (wxc.has("entropyFix"))
  {
    if (wxc.template get<std::string>("entropyFix") == "false")
      _efix = false;
  }
  // set minimum electron pressure to prevent negative
  if (wxc.has("minPressure"))
    _minPres = wxc.template get<REAL>("minPressure");
  else
    _minPres = 0.0;

  if (wxc.has("minDensity"))
    _minDens = wxc.template get<REAL>("minDensity");
  else
    _minDens = 0.0;

  // set Numerical flux to be used
  std::string lim;
  if (wxc.has("Numerical_Flux"))
      lim = wxc.template get<std::string>("Numerical_Flux");
  else
  {
    infoStrm << "WARNING: No Numerical Flux specified for Euler eqn., Lax-Friedrichs fluxes will be used.\n"
             << std::endl;

    lim = "LF";
    _fluxType = LF; // use LF flux evaluation
  }

  if (lim == "LF"){
    _fluxType = LF;
    infoStrm << "Using Lax-Friedrichs fluxes.\n"
             << std::endl;
  }
  else if (lim == "HLL"){
    _fluxType = HLL;
    infoStrm << "Using HLL fluxes.\n"
             << std::endl;
  }
  else if (lim == "Roe"){
    _fluxType = ROE;
    infoStrm << "Using Roe fluxes.\n"
             << std::endl;
  }
  else
  {
    infoStrm << "WARNING: Numerical FLux "
            << lim
            << " not recognised.\n Lax-Friedrichs fluxes will be used instead.\n"
            << std::endl;

    _fluxType = LF; // default to do Lax-Friedrichs fluxes.
  }

}

template<typename REAL>
void
WxEulerEntropyEqn<REAL>::
rp(unsigned d, REAL *ql, REAL *qr, REAL *qauxl, REAL *qauxr, REAL *df, REAL **wave, REAL *s, REAL *amdq, REAL *apdq)
{

}

template<typename REAL>
void
WxEulerEntropyEqn<REAL>::
rpt(unsigned td, unsigned d, REAL *ql, REAL* qr, REAL *amdq, REAL* bmamdq, REAL* bpamdq, REAL *apdq, REAL* bmapdq, REAL* bpapdq)
{

}

template<typename REAL>
void
WxEulerEntropyEqn<REAL>::
rptc(unsigned td, unsigned d, REAL *ql, REAL* qr,REAL *soc, REAL* bms, REAL* bps)
{

}

template<typename REAL>
void
WxEulerEntropyEqn<REAL>::
flux(unsigned d, REAL *x, REAL *q, REAL *qaux, REAL *f)
{
  double rho,u,v,w,S,p;
  unsigned mu=1, mv=2, mw=3;

  // depending on direction solve set velocity components. This is
  // done to ensure we are solving the Riemann problem across a face
  // in the direction 'd'
  if (d==0)
  { // x-direction Riemann problem
    mu = 1;
    mv = 2;
    mw = 3;
  }
  else if(d==1)
  { // y-direction Riemann problem
    mu = 2;
    mv = 3;
    mw = 1;
  }
  else if (d==2)
  { // z-direction Riemann problem
    mu = 3;
    mv = 1;
    mw = 2;
  }

  rho = q[0];
  if(_minDens == 0.0 && rho<=0){
      WxLogger *l = WxLogger::get("apollo-root.console");
      WxLogStream errStrm = l->getErrorStream();
      errStrm << "*** Negative density in Euler flux calculations. ***\n" ;
//      PetscFinalize();
      exit(1); // abort execution
  }
  if(rho < _minDens)
      rho = _minDens;

  u = q[mu]/rho;
  v = q[mv]/rho;
  w = q[mw]/rho;
  S = q[4];
  p = S*pow(rho,_gas_gamma);

  if(p<=0){
      if(_minPres==0.0){
          WxLogger *l = WxLogger::get("apollo-root.console");
          WxLogStream errStrm = l->getErrorStream();
          errStrm << "*** Negative pressure in Euler flux calculations. ***\n" ;
//          PetscFinalize();
          exit(1); // abort execution
      }
      else
          p = _minPres;
  }

  f[0] = rho*u;
  f[mu] = rho*u*u + p;
  f[mv] = rho*u*v;
  f[mw] = rho*u*w;
  f[4] = S*u;
}

template<typename REAL>
void
WxEulerEntropyEqn<REAL>::
DGnumericalFlux(REAL *normals, REAL *qM, REAL *qP, REAL *nflux, REAL *maxSpeed)
{
    switch (_fluxType) {
    case 0:
        applyLax_FriedrichsFluxes(normals,qM,qP,nflux,maxSpeed);
        break;
    case 1:
        applyHLLFluxes(normals,qM,qP,nflux,maxSpeed);
        break;
    case 2:
        applyRoeFluxes(normals,qM,qP,nflux,maxSpeed);
        break;
    default:
        applyLax_FriedrichsFluxes(normals,qM,qP,nflux,maxSpeed);
        break;
    }
}

template<typename REAL>
void
WxEulerEntropyEqn<REAL>::
fluxJacobian(unsigned d, REAL *x, REAL *q, REAL *qaux, REAL **f)
{

}

template<typename REAL>
void
WxEulerEntropyEqn<REAL>::
eigenSystem(unsigned d, REAL *q, REAL *ev, REAL **lev, REAL **rev)
{

}

template<typename REAL>
void
WxEulerEntropyEqn<REAL>::
RHS(unsigned N, REAL *geometry, REAL *normals, WxpDGGeometry<REAL> *quad, REAL *q, REAL *dq, REAL *rhs)
{

}

template<typename REAL>
void
WxEulerEntropyEqn<REAL>::
applyLax_FriedrichsFluxes(REAL *normals, REAL *qM, REAL *qP, REAL *nflux, REAL *maxSpeed)
{
    REAL *xc, *qaux;
    REAL fM[5], fP[5], gM[5], gP[5]; // x/y Fluxes
    REAL pM[5], pP[5]; // primitive variables

    // evaluate fluxes
    this->flux(0, xc, qM, qaux, fM);
    this->flux(0, xc, qP, qaux, fP);
    this->flux(1, xc, qM, qaux, gM);
    this->flux(1, xc, qP, qaux, gP);

    // compute primitive variables
    this->primitiveVariables(qM,pM);
    this->primitiveVariables(qP,pP);

    // compute the fastest propagating wave speed
    REAL c0M = sqrt(pM[1]*pM[1]+pM[2]*pM[2]+pM[3]*pM[3]) + sqrt(_gas_gamma*pM[4]/pM[0]);
    REAL c0P = sqrt(pP[1]*pP[1]+pP[2]*pP[2]+pP[3]*pP[3]) + sqrt(_gas_gamma*pP[4]/pP[0]);
    REAL lambda = dmax(c0M,c0P);

    // Lax-Frederick fluxes
    for(unsigned comp=0; comp<meqn(); comp++)
        nflux[comp] = 0.5*(normals[0]*(fM[comp]+fP[comp]) + normals[1]*(gM[comp]+gP[comp]) + lambda*(qM[comp]-qP[comp]));

    *maxSpeed = lambda;
}

template<typename REAL>
void
WxEulerEntropyEqn<REAL>::
applyHLLFluxes(REAL *normals, REAL *qM, REAL *qP, REAL *nflux, REAL *maxSpeed)
{
    REAL *xc, *qaux;
    REAL fM[5], fP[5], gM[5], gP[5]; // x/y Fluxes
    REAL pM[5], pP[5]; // primitive variables
    REAL fx[5];

//    if ((qM[0]<=0.) || (qP[0]<=0.))
//    {
//      WxLogger::get("apollo-root.console")->
//        error("*** Negative density in Euler HLL numerical flux calculation. *** ");
//      exit(1); // abort execution
//    }

    // Rotate "-" trace momentum to face normal-tangent coordinates
    REAL rhouM = qM[1], rhovM = qM[2];
    qM[1] = normals[0]*rhouM + normals[1]*rhovM;
    qM[2] =-normals[1]*rhouM + normals[0]*rhovM;

    // Rotate "+" trace momentum to face normal-tangent coordinates
    REAL rhouP = qP[1], rhovP = qP[2];
    qP[1] = normals[0]*rhouP + normals[1]*rhovP;
    qP[2] =-normals[1]*rhouP + normals[0]*rhovP;

    // evaluate fluxes and primitive variables in rotated coordinates
    this->flux(0, xc, qM, qaux, fM);
    this->flux(0, xc, qP, qaux, fP);
    this->flux(1, xc, qM, qaux, gM);
    this->flux(1, xc, qP, qaux, gP);

    // primitives
    this->primitiveVariables(qM,pM);
    this->primitiveVariables(qP,pP);

//    if ((pM[4]<=0.) || (pP[4]<=0.))
//    {
//      WxLogger::get("apollo-root.console")->
//        error("*** Negative pressure in Euler HLL numerical flux calculation. *** ");
//      exit(1); // abort execution
//    }

    REAL eM = qM[4]*pow(qM[0],_gas_gamma-1.)/(_gas_gamma-1.)+0.5*pM[0]*(pM[1]*pM[1]+pM[2]*pM[2]+pM[3]+pM[3]);
    REAL eP = qP[4]*pow(qP[0],_gas_gamma-1.)/(_gas_gamma-1.)+0.5*pP[0]*(pP[1]*pP[1]+pP[2]*pP[2]+pP[3]+pP[3]);

    REAL HM = (eM+pM[4])/pM[0], cM = sqrt(_gas_gamma*pM[4]/pM[0]);
    REAL HP = (eP+pP[4])/pP[0], cP = sqrt(_gas_gamma*pP[4]/pP[0]);

    // Compute Roe average variables
    REAL rhoMs = sqrt(pM[0]), rhoPs = sqrt(pP[0]);

    REAL u   = (rhoMs*pM[1] + rhoPs*pP[1])/(rhoMs + rhoPs);
    REAL v   = (rhoMs*pM[2] + rhoPs*pP[2])/(rhoMs + rhoPs);
    REAL w   = (rhoMs*pM[3] + rhoPs*pP[3])/(rhoMs + rhoPs);
    REAL H   = (rhoMs*HM    + rhoPs*HP)   /(rhoMs + rhoPs);

    REAL c2  = (_gas_gamma-1.)*(H - 0.5*(u*u + v*v +w*w)), c = sqrt(c2);

    // Compute estimate of waves speeds
    REAL SL = dmin(pM[1]-cM, u-c), SR = dmax(pP[1]+cP, u+c);

    // Compute HLL flux
    REAL zero = 0.0;
    REAL t1 = (dmin(SR,zero)-dmin(zero,SL))/(SR-SL);
    REAL t2 = 1.-t1;
    REAL t3 = (SR*fabs(SL)-SL*abs(SR))/(2.*(SR-SL));

    for(unsigned n=0; n<5; n++)
        fx[n] = t1*fP[n] + t2*fM[n] - t3*(qP[n]-qM[n]);

    // rotate flux back into Cartesian coordinates
    nflux[0] = fx[0];
    nflux[1] = normals[0]*fx[1] - normals[1]*fx[2];
    nflux[2] = normals[1]*fx[1] + normals[0]*fx[2];
    nflux[3] = 0.0;
    nflux[4] = fx[4];

    // compute the fastest propagating wave speed
    REAL c0M = sqrt(pM[1]*pM[1]+pM[2]*pM[2]+pM[3]*pM[3]) + sqrt(_gas_gamma*pM[4]/pM[0]);
    REAL c0P = sqrt(pP[1]*pP[1]+pP[2]*pP[2]+pP[3]*pP[3]) + sqrt(_gas_gamma*pP[4]/pP[0]);
    REAL lambda = dmax(c0M,c0P);

    *maxSpeed = lambda;
}

template<typename REAL>
void
WxEulerEntropyEqn<REAL>::
applyRoeFluxes(REAL *normals, REAL *qM, REAL *qP, REAL *nflux, REAL *maxSpeed)
{
    REAL *xc, *qaux;
    REAL fM[5], fP[5]; // x/y Fluxes
    REAL pM[5], pP[5]; // primitive variables
    REAL fx[5];

    // Rotate "-" trace momentum to face normal-tangent coordinates
    REAL rhouM = qM[1], rhovM = qM[2];
    qM[1] = normals[0]*rhouM + normals[1]*rhovM;
    qM[2] =-normals[1]*rhouM + normals[0]*rhovM;

    // Rotate "+" trace momentum to face normal-tangent coordinates
    REAL rhouP = qP[1], rhovP = qP[2];
    qP[1] = normals[0]*rhouP + normals[1]*rhovP;
    qP[2] =-normals[1]*rhouP + normals[0]*rhovP;

    // primitives
    this->primitiveVariables(qM,pM);
    this->primitiveVariables(qP,pP);

    REAL eM = qM[4]*pow(qM[0],_gas_gamma-1.)/(_gas_gamma-1.)+0.5*pM[0]*(pM[1]*pM[1]+pM[2]*pM[2]+pM[3]+pM[3]);
    REAL eP = qP[4]*pow(qP[0],_gas_gamma-1.)/(_gas_gamma-1.)+0.5*pP[0]*(pP[1]*pP[1]+pP[2]*pP[2]+pP[3]+pP[3]);

    REAL HM = (eM+pM[4])/pM[0];
    REAL HP = (eP+pP[4])/pP[0];

    // Compute Roe average variables
    REAL rhoMs = sqrt(pM[0]), rhoPs = sqrt(pP[0]);

    REAL rho = rhoMs*rhoPs;
    REAL u   = (rhoMs*pM[1] + rhoPs*pP[1])/(rhoMs + rhoPs);
    REAL v   = (rhoMs*pM[2] + rhoPs*pP[2])/(rhoMs + rhoPs);
    REAL w   = (rhoMs*pM[3] + rhoPs*pP[3])/(rhoMs + rhoPs);
    REAL H   = (rhoMs*HM    + rhoPs*HP)   /(rhoMs + rhoPs);

    REAL c2  = (_gas_gamma-1.)*(H - 0.5*(u*u + v*v +w*w)), c = sqrt(c2);

    // Riemann fluxes
    REAL dw1 = -0.5*rho*(pP[1]-pM[1])/c + 0.5*(pP[4]-pM[4])/c2;
    REAL dw2 = (pP[0]-pM[0]) - (pP[4]-pM[4])/c2;
    REAL dw3 = rho*(pP[2]-pM[2]);
    REAL dw4 = 0.5*rho*(pP[1]-pM[1])/c + 0.5*(pP[4]-pM[4])/c2;

    dw1 = fabs(u-c)*dw1;
    dw2 = fabs(u)*dw2;
    dw3 = fabs(u)*dw3;
    dw4 = fabs(u+c)*dw4;

    // From Roe fluxes
    // evaluate fluxes in rotated coordinates
    this->flux(0, xc, qM, qaux, fM);
    this->flux(0, xc, qP, qaux, fP);
    for(unsigned meqn=0; meqn<5; meqn++)
        fx[meqn] = 0.5*(fM[meqn]+fP[meqn]);

    fx[0] -= 0.5*(dw1*1.    + dw2*1. + dw3*0. + dw4*1.);
    fx[1] -= 0.5*(dw1*(u-c) + dw2*u  + dw3*0. + dw4*(u+c));
    fx[2] -= 0.5*(dw1*v     + dw2*v  + dw3*1. + dw4*v);
    fx[3]  = 0.0;
    fx[4] -= 0.5*(dw1*(H-u*c) + dw2*(u*u+v*v)/2. + dw3*v + dw4*(H+u*c));

    nflux[0] = fx[0];
    nflux[1] = normals[0]*fx[1] - normals[1]*fx[2];
    nflux[1] = normals[1]*fx[1] + normals[0]*fx[2];
    nflux[3] = 0.0;
    nflux[4] = fx[4];

    // compute the fastest propagating wave speed
    REAL c0M = sqrt(pM[1]*pM[1]+pM[2]*pM[2]+pM[3]*pM[3]) + sqrt(_gas_gamma*pM[4]/pM[0]);
    REAL c0P = sqrt(pP[1]*pP[1]+pP[2]*pP[2]+pP[3]*pP[3]) + sqrt(_gas_gamma*pP[4]/pP[0]);
    REAL lambda = dmax(c0M,c0P);

    *maxSpeed = lambda;
}

template<typename REAL>
void
WxEulerEntropyEqn<REAL>::
applyWavePropagationFluxes(REAL *normals, REAL *qM, REAL *qP, REAL *nflux, REAL *maxSpeed)
{

}

template<typename REAL>
void
WxEulerEntropyEqn<REAL>::
primitiveVariables(REAL *qCons, REAL *qPrim)
{
    qPrim[0] = qCons[0];
    if(qPrim[0]<=0.){
        qCons[0] = _minDens;
        qPrim[0] = _minDens;}

    qPrim[1] = qCons[1]/qCons[0];
    qPrim[2] = qCons[2]/qCons[0];
    qPrim[3] = qCons[3]/qCons[0];
    qPrim[4] = qCons[4]*pow(qCons[0],_gas_gamma);
}


template<typename REAL>
void
WxEulerEntropyEqn<REAL>::
limiterTuAndAliabadi(REAL *avgCons, REAL *avgPrim, REAL *dGrads, REAL *limitedValues)
{
    REAL Lrho, Lrhou, Lrhov, Lrhow, LEner, Lp;

    REAL averho = avgPrim[0];
    REAL aveu = avgPrim[1];
    REAL avev = avgPrim[2];
    REAL avew = avgPrim[3];
    REAL avep = avgPrim[4];

    REAL drho = dGrads[0];
    REAL du = dGrads[1];
    REAL dv = dGrads[2];
    REAL dw = dGrads[3];
    REAL dp = dGrads[4];

    Lrho = averho + drho;
    if(Lrho<_minDens)
    {
        for(unsigned riter=0; riter<4; riter++)
        {
            drho *= 0.5;
            Lrho = averho + drho;
            if (Lrho>_minDens)
                break;
            else if(riter==3)
              Lrho = avgCons[0];
        }
    }

    // Reconstruct momentum
    Lrhou = avgCons[1] + averho*du + drho*aveu;
    Lrhov = avgCons[2] + averho*dv + drho*avev;
    Lrhow = avgCons[3] + averho*dw + drho*avew;

    // Reconstruct energy
    REAL dEner = dp/(_gas_gamma-1.) + 0.5*drho*(aveu*aveu+avev*avev+avew*avew)
            + averho*(aveu*du+avev*dv+avew*dw);
    LEner = avgCons[4] + dEner;
    // Check for negative pressures => zero gradient
    Lp = (_gas_gamma-1.)*(LEner - 0.5*(Lrhou*Lrhou + Lrhov*Lrhov + Lrhow*Lrhow)/Lrho);

    if(Lp<_minPres)
    {
        for(unsigned piter=0; piter<4; piter++)
        {
            dp *= 0.5;
            dEner = dp/(_gas_gamma-1.) + 0.5*drho*(aveu*aveu+avev*avev+avew*avew)
                                + averho*(aveu*du+avev*dv+avew*dw);
            LEner = avgCons[4] + dEner;
            Lp = (_gas_gamma-1.)*(LEner - 0.5*(Lrhou*Lrhou + Lrhov*Lrhov + Lrhow*Lrhow)/Lrho);
            if (Lp>_minPres)
                break;
            else if(piter==3)
            {
                LEner = avep/(_gas_gamma-1.) + 0.5*(Lrhou*Lrhou+Lrhov*Lrhov+Lrhow*Lrhow)/Lrho;
            }
        }
    }

    limitedValues[0] = Lrho;
    if(Lrho!=Lrho){
        WxLogger::get("apollo-root.console")->
          error("*** NaN density in Euler limiter ***\n");
        exit(1); // abort execution
    }

    if(Lrho<=0.){
        WxLogger::get("apollo-root.console")->
          error("*** Negative density in Euler limiter ***\n");
        exit(1); // abort execution
    }

    limitedValues[1] = Lrhou;
    if(Lrhou!=Lrhou){
        WxLogger::get("apollo-root.console")->
          error("*** NaN u-velocity in Euler limiter ***\n");
        exit(1); // abort execution
    }

    limitedValues[2] = Lrhov;
    if(Lrhov!=Lrhov){
        WxLogger::get("apollo-root.console")->
          error("*** NaN v-velocity in Euler limiter ***\n");
        exit(1); // abort execution
    }

    limitedValues[3] = Lrhow;
    if(Lrho!=Lrho){
        WxLogger::get("apollo-root.console")->
          error("*** NaN w-velocity in Euler limiter ***\n");
        exit(1); // abort execution
    }

    limitedValues[4] = LEner;
    if(LEner!=LEner){
        WxLogger::get("apollo-root.console")->
          error("*** NaN energy in Euler limiter ***\n");
        exit(1); // abort execution
    }

//    if(Lp<=0.){
//        WxLogger::get("apollo-root.console")->
//          error("*** Negative pressure in Euler limiter ***\n");
//        exit(1); // abort execution
//    }
}

// instantiations
//template class WxEulerEntropyEqn<float>;
template class WxEulerEntropyEqn<double>;
