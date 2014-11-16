// WarpX lib includes
#include <wxlogger.h>

// WarpX hyperbolic solver includes
#include "wxeulereqn.h"
#include <wxlogger.h>
#include <wxlogstream.h>
// std includes
#include <cmath>

template<typename REAL>
void
WxEulerEqn<REAL>::
setup(const WxCryptSet& wxc)
{
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

}

template<typename REAL>
void
WxEulerEqn<REAL>::
rp(unsigned d, REAL *ql, REAL *qr, REAL *qauxl, REAL *qauxr, REAL *df, REAL **wave, REAL *s, REAL *amdq, REAL *apdq)
{
  // This solver is based on Roe-averages. The waves with same
  // eigenvalues are lumped into a single wave.

  REAL rhsqrtl,rhsqrtr,pl,pr,rhsq2;
  REAL a1,a2,a3,a4,a5;
  REAL u,v,w,u2v2w2,enth,aa2,a,g1a2,euv;
  REAL gas_gamma = _gas_gamma;
  REAL gas_gamma1 = _gas_gamma-1;
  REAL delta[5];
  int* idx;

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

  //
  // Compute Roe-averaged quantities
  //

  if ((ql[0]<0) || (qr[0]<0))
  {
//		this -> getIndices(idx);
		std::stringstream ss;
		WxLogger *l = WxLogger::get("warpx-root.console");
		WxLogStream errStrm = l->getErrorStream();
//		errStrm << "*** Negative or zero density in Euler Riemann solver at index = (" << idx[0] << "," << idx[1] <<  "," <<  idx[2] <<  ")" ;
		errStrm << "*** Negative or zero density in Euler Riemann solver" ;
	    exit(1); // abort execution
  }
  rhsqrtl = sqrt(ql[0]);
  rhsqrtr = sqrt(qr[0]);
  // left edge pressure
  pl = gas_gamma1*(ql[4]
    - 0.5*(dsqr(ql[mu]) + dsqr(ql[mv]) + dsqr(ql[mw]))/ql[0]);
  // right edge pressure
  pr = gas_gamma1*(qr[4]
    - 0.5*(dsqr(qr[mu]) + dsqr(qr[mv]) + dsqr(qr[mw]))/qr[0]);
  if (_minPres==0.0 && (pl<0 || pr<0))
  {
//		this -> getIndices(idx);
		std::stringstream ss;
		WxLogger *l = WxLogger::get("warpx-root.console");
		WxLogStream errStrm = l->getErrorStream();
//		errStrm << "*** Negative or zero pressure in Euler Riemann solver at index = (" << idx[0] << "," << idx[1] <<  "," <<  idx[2] <<  ")" ;
		errStrm << "*** Negative or zero pressure in Euler Riemann solver" ;
		exit(1); // abort execution
  }    

  // if min pressures specified then set the floor values
  if (pl<_minPres) 
  {
    ql[4] = _minPres/(_gas_gamma-1) + 0.5*(dsqr(ql[1]) + dsqr(ql[2]) + dsqr(ql[3]))/ql[0];
    pl = _minPres;
  }
  if (pr<_minPres) 
  {
    qr[4] = _minPres/(_gas_gamma-1) + 0.5*(dsqr(qr[1]) + dsqr(qr[2]) + dsqr(qr[3]))/qr[0];
    pr = _minPres;
  }

  rhsq2 = rhsqrtl + rhsqrtr;
  // Roe-averaged velocity components
  u = (ql[mu]/rhsqrtl + qr[mu]/rhsqrtr) / rhsq2;
  v = (ql[mv]/rhsqrtl + qr[mv]/rhsqrtr) / rhsq2;
  w = (ql[mw]/rhsqrtl + qr[mw]/rhsqrtr) / rhsq2;
  u2v2w2 = u*u + v*v + w*w;
  // Roe-averaged enthalpy
  enth = ( (ql[4]+pl)/rhsqrtl + (qr[4]+pr)/rhsqrtr ) / rhsq2;
  // speed of sound
  aa2 = gas_gamma1*(enth - 0.5*u2v2w2);
  a = sqrt(aa2);
  if(aa2<0)
  {
//		this -> getIndices(idx);
		std::stringstream ss;
		WxLogger *l = WxLogger::get("warpx-root.console");
		WxLogStream errStrm = l->getErrorStream();
//		errStrm << "*** Negative sound-speed in Euler Riemann solver at index = (" << idx[0] << "," << idx[1] <<  "," <<  idx[2] <<  ")" ;
		errStrm << "*** Negative sound-speed in Euler Riemann solver at index" ;
		exit(1); // abort execution

  }
  g1a2 = gas_gamma1 / aa2;
  euv  = enth - u2v2w2;

  //
  // Compute waves
  //

  delta[0] = df[0];
  delta[1] = df[mu];
  delta[2] = df[mv];
  delta[3] = df[mw];
  delta[4] = df[4];

  // calculate coefficients of the 5 eigenvectors
  a4 = g1a2 * (euv*delta[0] + u*delta[1] + v*delta[2] + w*delta[3] - delta[4]);
  a2 = delta[2] - v*delta[0];
  a3 = delta[3] - w*delta[0];
  a5 = (delta[1] + (a-u)*delta[0] - a*a4) / (2.0*a);
  a1 = delta[0] - a4 - a5;
    
  // Wave 1: eigenvalue is u-c
  wave[0][0]  = a1;
  wave[mu][0] = a1*(u-a);
  wave[mv][0] = a1*v;
  wave[mw][0] = a1*w;
  wave[4][0]  = a1*(enth - u*a);
  s[0] = u-a;
    
  // Wave 2: the 3 eigenvectors corresponding to the repeated
  // eigenvalue u,u,u are lumped together into a single wave
  wave[0][1]  = a4;
  wave[mu][1] = a4*u;
  wave[mv][1] = a4*v	 	 + a2;
  wave[mw][1] = a4*w	 	 + a3;
  wave[4][1]  = a4*0.5*u2v2w2  + a2*v + a3*w;
  s[1] = u;
    
  // Wave 3: eigenvalue is u+c
  wave[0][2]  = a5;
  wave[mu][2] = a5*(u+a);
  wave[mv][2] = a5*v;
  wave[mw][2] = a5*w;
  wave[4][2]  = a5*(enth+u*a);
  s[2] = u+a;

  if (!_efix) 
  {
    // compute fluctuations
    for (unsigned m=0; m<5; ++m)
    {
      amdq[m] = 0.0; apdq[m] = 0.0;
      for (unsigned mw=0; mw<3; ++mw)
      {
        if (s[mw] < 0.0)
          // left going wave
          amdq[m] += s[mw]*wave[m][mw];
        else
          // right going wave
          apdq[m] += s[mw]*wave[m][mw];
      }
    }
    // return right away if we are not applying entropy fix
    return;
  }

  REAL rhoim1, pim1, cim1, s0;
  REAL rho1, rhou1, rhov1, rhow1, en1, p1, c1, s1;
  REAL rhoi, pi, ci, s3;
  REAL rho2, rhou2, rhov2, rhow2, en2, p2, c2, s2;
  REAL sfract, ddf;

  // apply entropy fix
  rhoim1 = ql[0];
  pim1 = gas_gamma1*(ql[4] - 0.5*(ql[mu]*ql[mu] + ql[mv]*ql[mv] + ql[mw]*ql[mw])/rhoim1);
  cim1 = sqrt(gas_gamma*pim1/rhoim1);
  s0 = ql[mu]/rhoim1 - cim1; // u-c in left state (cell i-1)

  // check for fully supersonic case:
  if ((s0 > 0) &&  s[0] > 0)
  {
    // everything is right-going
    for (unsigned m=0; m<5; ++m)
      amdq[m] = 0.0;
    goto computeapdq;
  }

  rho1 = ql[0] + wave[0][0];
  rhou1 = ql[mu] + wave[mu][0];
  rhov1 = ql[mv] + wave[mv][0];
  rhow1 = ql[mw] + wave[mw][0];
  en1 = ql[4] + wave[4][0];
  p1 = gas_gamma1*(en1 - 0.5*(rhou1*rhou1 + rhov1*rhov1 + rhow1*rhow1)/rho1);
  c1 = sqrt(gas_gamma*p1/rho1);
  s1 = rhou1/rho1 - c1; // u-c to right of 1-wave

  if ((s0<0.0) && (s1>0.0))
    // transonic rarefaction in the 1-wave
    sfract = s0*(s1-s[0]) / (s1-s0);
  else if (s[0] < 0.0)
    // 1-wave is leftgoing
    sfract = s[0];
  else
    // 1-wave is rightgoing
    sfract = 0.0; // this shouldn't happen since s0 < 0

  for (unsigned m=0; m<5; ++m)
    amdq[m] = sfract*wave[m][0];

  // check 2-wave:

  if (s[1] >= 0.0) 
    goto computeapdq; // 2-,3- and 4- waves are rightgoing
  for (unsigned m=0; m<5; ++m)
    amdq[m] = amdq[m] + s[1]*wave[m][1];

  // check 3-wave:
  rhoi = qr[0];
  pi = gas_gamma1*(qr[4] - 0.5*(qr[mu]*qr[mu] + qr[mv]*qr[mv] + qr[mw]*qr[mw])/rhoi);
  ci = sqrt(gas_gamma*pi/rhoi);
  s3 = qr[mu]/rhoi + ci; // u+c in right state (cell i)

  rho2 = qr[0] - wave[0][2];
  rhou2 = qr[mu] - wave[mu][2];
  rhov2 = qr[mv] - wave[mv][2];
  rhow2 = qr[mw] - wave[mw][2];
  en2 = qr[4] - wave[4][2];
  p2 = gas_gamma1*(en2 - 0.5*(rhou2*rhou2 + rhov2*rhov2 + rhow2*rhow2)/rho2);
  c2 = sqrt(gas_gamma*p2/rho2);
  s2 = rhou2/rho2 + c2; // u+c to left of 3-wave
  if ((s2 < 0.0) && (s3 > 0.0))
    // transonic rarefaction in the 3-wave
    sfract = s2 * (s3-s[2]) / (s3-s2);
  else if (s[2] < 0.0)
    // 3-wave is leftgoing
    sfract = s[2];
  else
    // 3-wave is rightgoing
    goto computeapdq;

  for (unsigned m=0; m<5; ++m)
    amdq[m] = amdq[m] + sfract*wave[m][2];

  computeapdq:
  ;
    
  // compute the rightgoing flux differences:
  // ddf = SUM s*wave   is the total flux difference and apdq = ddf - amdq

  for (unsigned m=0; m<5; ++m)
  {
    ddf = 0.0;
    for (unsigned mw=0; mw<3; ++mw)
      ddf += s[mw]*wave[m][mw];
    apdq[m] = ddf - amdq[m];
  }
}

template<typename REAL>
void
WxEulerEqn<REAL>::
rpt(unsigned td, unsigned d, REAL *ql, REAL* qr, REAL *amdq, REAL* bmamdq, REAL* bpamdq, REAL *apdq, REAL* bmapdq, REAL* bpapdq)
{
  REAL rhsqrtl,rhsqrtr,pl,pr,rhsq2;
  REAL a1,a2,a3,a4,a5;
  REAL u,v,w,u2v2w2,enth,aa2,a,g1a2,euv;
  //REAL gas_gamma = _gas_gamma;
  REAL gas_gamma1 = _gas_gamma-1;
  REAL delta[5];
  REAL wave[5][3];
  REAL s[3];
  REAL vel[3];
  unsigned mu=1,  mv=2,  mw=3;
  unsigned mmu=1, mmv=2, mmw=3;
  int* idx;


  // depending on direction solve set velocity components. This is
  // done to ensure we are solving the Riemann problem across a face
  // in the direction 'd'
  if (d==0)
  { // x-direction Riemann problem
	  mmu = 1;
	  mmv = 2;
	  mmw = 3;

	  if (td==1)
	  {
		  mu = 2;
	  	  mv = 3;
	  	  mw = 1;
	  }
	  else if (td==2)
	  {
		  mu = 3;
	  	  mv = 1;
	  	  mw = 2;
	  }
	  else
		  throw WxExcept("WxEulerEqn::rpt: transverse direction code not valid");

  }
  else if(d==1)
  { // y-direction Riemann problem
	  mmu = 2;
	  mmv = 3;
	  mmw = 1;

	  if (td==2)
	  {
		  mu = 3;
	  	  mv = 1;
	  	  mw = 2;
	  }
	  else if (td==0)
	  {
		  mu = 1;
	  	  mv = 2;
	  	  mw = 3;
	  }
	  else
		  throw WxExcept("WxEulerEqn::rpt: transverse direction code not valid");
  }

  else if (d==2)
  { // z-direction Riemann problem

    mmu = 3;
    mmv = 1;
    mmw = 2;

    if (td==0)
    {
    	mu = 1;
    	mv = 2;
    	mw = 3;
    }
    else if (td==1)
    {
    	mu = 2;
    	mv = 3;
    	mw = 1;
    }
    else
    	throw WxExcept("WxEulerEqn::rpt: transverse direction code not valid");
  }

  rhsqrtl = sqrt(ql[0]);
  rhsqrtr = sqrt(qr[0]);

  // left edge pressure
//  pl = gas_gamma1*(ql[4] - 0.5*(dsqr(ql[mu]) + dsqr(ql[mv]) + dsqr(ql[mw]))/ql[0]);
    pl = gas_gamma1*(ql[4] - 0.5*(dsqr(ql[mmu]) + dsqr(ql[mmv]) + dsqr(ql[mmw]))/ql[0]);
  // right edge pressure
//  pr = gas_gamma1*(qr[4] - 0.5*(dsqr(qr[mu]) + dsqr(qr[mv]) + dsqr(qr[mw]))/qr[0]);
    pr = gas_gamma1*(qr[4] - 0.5*(dsqr(qr[mmu]) + dsqr(qr[mmv]) + dsqr(qr[mmw]))/qr[0]);

  if ((pl<0) || (pr<0))
  {
//	  this -> getIndices(idx);
	  std::stringstream ss;
	  WxLogger *l = WxLogger::get("warpx-root.console");
	  WxLogStream errStrm = l->getErrorStream();
//	  errStrm << "*** Negative or zero pressure in Euler in Euler Transverse Riemann solver at index = (" << idx[0] << "," << idx[1] <<  "," <<  idx[2] <<  ")" ;
	  errStrm << "*** Negative or zero pressure in Euler in Euler Transverse Riemann solver" ;
	  exit(1); // abort execution
  }

  rhsq2 = rhsqrtl + rhsqrtr;

  // Roe-averaged velocity components
    vel[mmu] = (ql[mmu]/rhsqrtl + qr[mmu]/rhsqrtr) / rhsq2;
    vel[mmv] = (ql[mmv]/rhsqrtl + qr[mmv]/rhsqrtr) / rhsq2;
    vel[mmw] = (ql[mmw]/rhsqrtl + qr[mmw]/rhsqrtr) / rhsq2;

//  u = (ql[mu]/rhsqrtl + qr[mu]/rhsqrtr) / rhsq2;
//  v = (ql[mv]/rhsqrtl + qr[mv]/rhsqrtr) / rhsq2;
//  w = (ql[mw]/rhsqrtl + qr[mw]/rhsqrtr) / rhsq2;

    u = vel[mu];
    v = vel[mv];
    w = vel[mw];

  u2v2w2 = u*u + v*v + w*w;

  // Roe-averaged enthalpy
  enth = ( (ql[4]+pl)/rhsqrtl + (qr[4]+pr)/rhsqrtr ) / rhsq2;
  // speed of sound
  aa2 = gas_gamma1*(enth - 0.5*u2v2w2);
  a = sqrt(aa2);

  if(aa2<0)
  {
//		this -> getIndices(idx);
		std::stringstream ss;
		WxLogger *l = WxLogger::get("warpx-root.console");
		WxLogStream errStrm = l->getErrorStream();
//		errStrm << "*** Negative sound-speed in Euler Transverse Riemann solver at index = (" << idx[0] << "," << idx[1] <<  "," <<  idx[2] <<  ")" ;
		errStrm << "*** Negative sound-speed in Euler Transverse Riemann solver" ;
	    exit(1); // abort execution
  }

    g1a2 = gas_gamma1 / aa2;
  euv  = enth - u2v2w2;

  delta[0] = amdq[0];
  delta[1] = amdq[mu];
  delta[2] = amdq[mv];
  delta[3] = amdq[mw];
  delta[4] = amdq[4];

  // calculate coefficients of the 5 eigenvectors
  a4 = g1a2 * (euv*delta[0] + u*delta[1] + v*delta[2] + w*delta[3] - delta[4]);
  a2 = delta[2] - v*delta[0];
  a3 = delta[3] - w*delta[0];
  a5 = (delta[1] + (a-u)*delta[0] - a*a4) / (2.0*a);
  a1 = delta[0] - a4 - a5;

  // Wave 1: eigenvalue is u-c
  wave[0][0]  = a1;
  wave[mu][0] = a1*(u-a);
  wave[mv][0] = a1*v;
  wave[mw][0] = a1*w;
  wave[4][0]  = a1*(enth - u*a);
  s[0] = u-a;

  // Wave 2: the 3 eigenvectors corresponding to the repeated
  // eigenvalue u,u,u are lumped together into a single wave
  wave[0][1]  = a4;
  wave[mu][1] = a4*u;
  wave[mv][1] = a4*v	       + a2;
  wave[mw][1] = a4*w	 	      + a3;
  wave[4][1]  = a4*0.5*u2v2w2  + a2*v + a3*w;
  s[1] = u;

  // Wave 3: eigenvalue is u+c
  wave[0][2]  = a5;
  wave[mu][2] = a5*(u+a);
  wave[mv][2] = a5*v;
  wave[mw][2] = a5*w;
  wave[4][2]  = a5*(enth+u*a);
  s[2] = u+a;

  for (unsigned m=0; m<5; ++m)
  {
    bmamdq[m] = 0.0; bpamdq[m] = 0.0;
    for (unsigned mw=0; mw<3; ++mw)
    {
      if (s[mw] < 0.0)
        // left going wave
        bmamdq[m] += s[mw]*wave[m][mw];
      else
        // right going wave
        bpamdq[m] += s[mw]*wave[m][mw];
    }
  }

  delta[0] = apdq[0];
  delta[1] = apdq[mu];
  delta[2] = apdq[mv];
  delta[3] = apdq[mw];
  delta[4] = apdq[4];

  // calculate coefficients of the 5 eigenvectors
  //a4 = g1a2 * (euv*delta[0] + u*delta[1] + v*delta[2] + w*delta[3] - delta[4]);
  //a2 = delta[1] - u*delta[0];
  //a3 = delta[3] - w*delta[0];
  //a5 = (delta[2] + (a-v)*delta[0] - a*a4) / (2.0*a);
  //a1 = delta[0] - a4 - a5;

  // calculate coefficients of the 5 eigenvectors
  a4 = g1a2 * (euv*delta[0] + u*delta[1] + v*delta[2] + w*delta[3] - delta[4]);
  a2 = delta[2] - v*delta[0];
  a3 = delta[3] - w*delta[0];
  a5 = (delta[1] + (a-u)*delta[0] - a*a4) / (2.0*a);
  a1 = delta[0] - a4 - a5;

  // Wave 1: eigenvalue is u-c
  wave[0][0]  = a1;
  wave[mu][0] = a1*(u-a);
  wave[mv][0] = a1*v;
  wave[mw][0] = a1*w;
  wave[4][0]  = a1*(enth - u*a);
  s[0] = u-a;

  // Wave 2: the 3 eigenvectors corresponding to the repeated
  // eigenvalue u,u,u are lumped together into a single wave
  wave[0][1]  = a4;
  wave[mu][1] = a4*u;
  wave[mv][1] = a4*v	 	 + a2;
  wave[mw][1] = a4*w	 	 + a3;
  wave[4][1]  = a4*0.5*u2v2w2  + a2*v + a3*w;
  s[1] = u;

  // Wave 3: eigenvalue is u+c
  wave[0][2]  = a5;
  wave[mu][2] = a5*(u+a);
  wave[mv][2] = a5*v;
  wave[mw][2] = a5*w;
  wave[4][2]  = a5*(enth+u*a);
  s[2] = u+a;

  for (unsigned m=0; m<5; ++m)
  {
    bmapdq[m] = 0.0; bpapdq[m] = 0.0;
    for (unsigned mw=0; mw<3; ++mw)
    {
      if (s[mw] < 0.0)
        // left going wave
        bmapdq[m] += s[mw]*wave[m][mw];
      else
        // right going wave
        bpapdq[m] += s[mw]*wave[m][mw];
    }
  }
}

template<typename REAL>
void
WxEulerEqn<REAL>::
rptc(unsigned td, unsigned d, REAL *ql, REAL* qr,REAL *soc, REAL* bms, REAL* bps)
{
  REAL rhsqrtl,rhsqrtr,pl,pr,rhsq2;
  REAL a1,a2,a3,a4,a5;
  REAL u,v,w,u2v2w2,enth,aa2,a,g1a2,euv;
  //REAL gas_gamma = _gas_gamma;
  REAL gas_gamma1 = _gas_gamma-1;
  REAL delta[5];
  REAL wave[5][3];
  REAL s[3];
  REAL vel[3];
  unsigned mu=1,  mv=2,  mw=3;
  unsigned mmu=1, mmv=2, mmw=3;
  int* idx;

  // depending on direction solve set velocity components. This is
  // done to ensure we are solving the Riemann problem across a face
  // in the direction 'd'
  if (d==0)
  { // x-direction Riemann problem
	  mmu = 1;
	  mmv = 2;
	  mmw = 3;

	  if (td==1)
	  {
		  mu = 2;
	  	  mv = 3;
	  	  mw = 1;
	  }
	  else if (td==2)
	  {
		  mu = 3;
	  	  mv = 1;
	  	  mw = 2;
	  }
	  else
		  throw WxExcept("WxEulerEqn::rpt: transverse direction code not valid");

  }
  else if(d==1)
  { // y-direction Riemann problem
	  mmu = 2;
	  mmv = 3;
	  mmw = 1;

	  if (td==2)
	  {
		  mu = 3;
	  	  mv = 1;
	  	  mw = 2;
	  }
	  else if (td==0)
	  {
		  mu = 1;
	  	  mv = 2;
	  	  mw = 3;
	  }
	  else
		  throw WxExcept("WxEulerEqn::rpt: transverse direction code not valid");
  }

  else if (d==2)
  { // z-direction Riemann problem

    mmu = 3;
    mmv = 1;
    mmw = 2;

    if (td==0)
    {
    	mu = 1;
    	mv = 2;
    	mw = 3;
    }
    else if (td==1)
    {
    	mu = 2;
    	mv = 3;
    	mw = 1;
    }
    else
    	throw WxExcept("WxEulerEqn::rpt: transverse direction code not valid");
  }

  rhsqrtl = sqrt(ql[0]);
  rhsqrtr = sqrt(qr[0]);

  // left edge pressure
//  pl = gas_gamma1*(ql[4] - 0.5*(dsqr(ql[mu]) + dsqr(ql[mv]) + dsqr(ql[mw]))/ql[0]);
    pl = gas_gamma1*(ql[4] - 0.5*(dsqr(ql[mmu]) + dsqr(ql[mmv]) + dsqr(ql[mmw]))/ql[0]);
  // right edge pressure
//  pr = gas_gamma1*(qr[4] - 0.5*(dsqr(qr[mu]) + dsqr(qr[mv]) + dsqr(qr[mw]))/qr[0]);
    pr = gas_gamma1*(qr[4] - 0.5*(dsqr(qr[mmu]) + dsqr(qr[mmv]) + dsqr(qr[mmw]))/qr[0]);

  if ((pl<0) || (pr<0))
  {
//		this -> getIndices(idx);
		std::stringstream ss;
		WxLogger *l = WxLogger::get("warpx-root.console");
		WxLogStream errStrm = l->getErrorStream();
//		errStrm << "*** Negative pressure in Euler Transverse Correction Wave Riemann solver at index = (" << idx[0] << "," << idx[1] <<  "," <<  idx[2] <<  ")" ;
		errStrm << "*** Negative pressure in Euler Transverse Correction Wave Riemann solver at index" ;
	    exit(1); // abort execution
  }

  rhsq2 = rhsqrtl + rhsqrtr;

  // Roe-averaged velocity components
    vel[mmu] = (ql[mmu]/rhsqrtl + qr[mmu]/rhsqrtr) / rhsq2;
    vel[mmv] = (ql[mmv]/rhsqrtl + qr[mmv]/rhsqrtr) / rhsq2;
    vel[mmw] = (ql[mmw]/rhsqrtl + qr[mmw]/rhsqrtr) / rhsq2;

    u = vel[mu];
    v = vel[mv];
    w = vel[mw];

//  u = (ql[mu]/rhsqrtl + qr[mu]/rhsqrtr) / rhsq2;
//  v = (ql[mv]/rhsqrtl + qr[mv]/rhsqrtr) / rhsq2;
//  w = (ql[mw]/rhsqrtl + qr[mw]/rhsqrtr) / rhsq2;


  u2v2w2 = u*u + v*v + w*w;

  // Roe-averaged enthalpy
  enth = ( (ql[4]+pl)/rhsqrtl + (qr[4]+pr)/rhsqrtr ) / rhsq2;
  // speed of sound
  aa2 = gas_gamma1*(enth - 0.5*u2v2w2);
  a = sqrt(aa2);

  if(aa2<0)
  {
//	  this -> getIndices(idx);
	  std::stringstream ss;
	  WxLogger *l = WxLogger::get("warpx-root.console");
	  WxLogStream errStrm = l->getErrorStream();
//	  errStrm << "*** Negative sound-speed in Euler Transverse Correction Wave Riemann solver at index = (" << idx[0] << "," << idx[1] <<  "," <<  idx[2] <<  ")" ;
	  errStrm << "*** Negative sound-speed in Euler Transverse Correction Wave Riemann solver" ;
	  exit(1); // abort execution
  }

    g1a2 = gas_gamma1 / aa2;
  euv  = enth - u2v2w2;

  delta[0] = soc[0];
  delta[1] = soc[mu];
  delta[2] = soc[mv];
  delta[3] = soc[mw];
  delta[4] = soc[4];

  // calculate coefficients of the 5 eigenvectors
  a4 = g1a2 * (euv*delta[0] + u*delta[1] + v*delta[2] + w*delta[3] - delta[4]);
  a2 = delta[2] - v*delta[0];
  a3 = delta[3] - w*delta[0];
  a5 = (delta[1] + (a-u)*delta[0] - a*a4) / (2.0*a);
  a1 = delta[0] - a4 - a5;

  // Wave 1: eigenvalue is u-c
  wave[0][0]  = a1;
  wave[mu][0] = a1*(u-a);
  wave[mv][0] = a1*v;
  wave[mw][0] = a1*w;
  wave[4][0]  = a1*(enth - u*a);
  s[0] = u-a;

  // Wave 2: the 3 eigenvectors corresponding to the repeated
  // eigenvalue u,u,u are lumped together into a single wave
  wave[0][1]  = a4;
  wave[mu][1] = a4*u;
  wave[mv][1] = a4*v	       + a2;
  wave[mw][1] = a4*w	 	      + a3;
  wave[4][1]  = a4*0.5*u2v2w2  + a2*v + a3*w;
  s[1] = u;

  // Wave 3: eigenvalue is u+c
  wave[0][2]  = a5;
  wave[mu][2] = a5*(u+a);
  wave[mv][2] = a5*v;
  wave[mw][2] = a5*w;
  wave[4][2]  = a5*(enth+u*a);
  s[2] = u+a;

  for (unsigned m=0; m<5; ++m)
  {
    bms[m] = 0.0; bps[m] = 0.0;
    for (unsigned mw=0; mw<3; ++mw)
    {
      if (s[mw] < 0.0)
        // left going wave
        bms[m] += s[mw]*wave[m][mw];
      else
        // right going wave
        bps[m] += s[mw]*wave[m][mw];
    }
  }
}

template<typename REAL>
void
WxEulerEqn<REAL>::
flux(unsigned d, REAL *x, REAL *q, REAL *qaux, REAL *f)
{
  double rho,u,v,w,E,p;
  double gas_gamma1 = _gas_gamma-1.;
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
  u = q[mu]/rho;
  v = q[mv]/rho;
  w = q[mw]/rho;
  E = q[4];
  p = gas_gamma1*(E-0.5*rho*(u*u+v*v+w*w));

  f[0] = rho*u;
  f[mu] = rho*u*u + p;
  f[mv] = rho*u*v;
  f[mw] = rho*u*w;
  f[4] = (E+p)*u;
}

template<typename REAL>
void
WxEulerEqn<REAL>::
fluxJacobian(unsigned d, REAL *x, REAL *q, REAL *qaux, REAL **f)
{
  double gm1 = _gas_gamma-1.;
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


  f[0][0] = 0.0;
  f[0][mu]= 1.0;
  f[0][mv]= 0.0;
  f[0][mw]= 0.0;
  f[0][4] = 0.0;

  f[mu][0] = -(q[mu]*q[mu])/(q[0]*q[0]) + gm1*(q[mu]*q[mu]+q[mv]*q[mv]+q[mw]*q[mw])/(2.*q[0]*q[0]);
  f[mu][mu]= 2*q[mu]/q[0] - gm1*q[mu]/q[0];
  f[mu][mv]= -gm1*q[mv]/q[0];
  f[mu][mw]= -gm1*q[mw]/q[0];
  f[mu][4] = gm1;

  f[mv][0] = -q[mu]*q[mv]/(q[0]*q[0]);
  f[mv][mu]= q[mv]/q[0];
  f[mv][mv]= q[mu]/q[0];
  f[mv][mw]= 0.0;
  f[mv][4] = 0.0;

  f[mw][0] = -q[mu]*q[mw]/(q[0]*q[0]);
  f[mw][mu]= q[mw]/q[0];
  f[mw][mv]= 0.0;
  f[mw][mw]= q[mu]/q[0];
  f[mw][4] = 0.0;

  f[4][0] = -q[mu]*q[4]/(q[0]*q[0])-q[mu]*gm1/(q[0]*q[0])*(q[4]-(q[mu]*q[mu]+q[mv]*q[mv]+q[mw]*q[mw])/(2*q[0]))+gm1*q[mu]/q[0]*(q[mu]*q[mu]+q[mv]*q[mv]+q[mw]*q[mw])/(2*q[0]*q[0]);
  f[4][mu]= q[4]/q[0]+gm1/q[0]*(q[4]-(q[mu]*q[mu]+q[mv]*q[mv]+q[mw]*q[mw])/(2*q[0]))-gm1*q[mu]*q[mu]/(q[0]*q[0]);
  f[4][mv]= -gm1*q[mu]*q[mv]/(q[0]*q[0]);
  f[4][mw]= -gm1*q[mu]*q[mw]/q[0]/q[0];
  f[4][4] = q[mu]/q[0]+gm1*q[mu]/q[0];
}

template<typename REAL>
void
WxEulerEqn<REAL>::
eigenSystem(unsigned d, REAL *q, REAL *ev, REAL **lev, REAL **rev)
{
  unsigned mu=1, mv=2, mw=3;
  if (d==0)
  { // x-direction eigensystem
    mu = 1;
    mv = 2;
    mw = 3;
  }
  else if(d==1)
  { // y-direction eigensystem
    mu = 2;
    mv = 3;
    mw = 1;
  }
  else if (d==2)
  { // z-direction eigensystem
    mu = 3;
    mv = 1;
    mw = 2;
  }

  double rho, u, v, w, E, t, p, gas_gamma1, H;
  gas_gamma1 = _gas_gamma - 1;

  // compute conserved variables
  rho = q[0];
  u = q[mu]/rho;
  v = q[mv]/rho;
  w = q[mw]/rho;
  E = q[4];
  t = u*u + v*v + w*w;
  p = gas_gamma1*(E-0.5*rho*t);

  if (p<_minPres && _minPres==0.0)
  {
    std::cout<<"p = "<<p;
    WxLogger::get("warpx-root.console")->
      error("*** Negative pressure in Euler eigenSystem");
    exit(1); // abort execution
  }    
  if (p<_minPres)
    p = _minPres;

  H = (E+p)/rho;

  // sound speed
  REAL c = sqrt(_gas_gamma*p/rho);
    
  // eigenvalues
  ev[0] = u-c;
  ev[1] = ev[2] = ev[3] = u;
  ev[4] = u+c;

  // right eigenvectors: these are arranged as column vectors
  rev[0][0]  = 1.0;
  rev[mu][0] = u-c;
  rev[mv][0] = v;
  rev[mw][0] = w;
  rev[4][0]  = H - u*c;

  rev[0][1]  = 0.0;
  rev[mu][1] = 0.0;
  rev[mv][1] = 1.0;
  rev[mw][1] = 0.0;
  rev[4][1]  = v;

  rev[0][2]  = 1.0;
  rev[mu][2] = u;
  rev[mv][2] = v;
  rev[mw][2] = w;
  rev[4][2]  = 0.5*t;

  rev[0][3]  = 0.0;
  rev[mu][3] = 0.0;
  rev[mv][3] = 0.0;
  rev[mw][3] = 1.0;
  rev[4][3]  = w;

  rev[0][4]  = 1.0;
  rev[mu][4] = u+c;
  rev[mv][4] = v;
  rev[mw][4] = w;
  rev[4][4]  = H + u*c;

  REAL c2 = 1/(c*c);
  REAL c2h = 0.5*c2;
  // left eigenvectors: these are arranged as row vectors
  lev[0][0]  = c2h*(0.5*gas_gamma1*t + c*u);
  lev[0][mu] = -c2h*(gas_gamma1*u + c);
  lev[0][mv] = -c2h*gas_gamma1*v;
  lev[0][mw] = -c2h*gas_gamma1*w;
  lev[0][4]  = c2h*gas_gamma1;

  lev[1][0]  = -v;
  lev[1][mu] = 0.0;
  lev[1][mv] = 1.0;
  lev[1][mw] = 0.0;
  lev[1][4]  = 0.0;

  lev[2][0]  = 2 - gas_gamma1*H*c2;
  lev[2][mu] = c2*gas_gamma1*u;
  lev[2][mv] = c2*gas_gamma1*v;
  lev[2][mw] = c2*gas_gamma1*w;
  lev[2][4]  = -c2*gas_gamma1;

  lev[3][0] = -w;
  lev[3][mu] = 0.0;
  lev[3][mv] = 0.0;
  lev[3][mw] = 1.0;
  lev[3][4] = 0.0;

  lev[4][0]  = c2h*(0.5*gas_gamma1*t - c*u);
  lev[4][mu] = -c2h*(gas_gamma1*u - c);
  lev[4][mv] = -c2h*gas_gamma1*v;
  lev[4][mw] = -c2h*gas_gamma1*w;
  lev[4][4]  = c2h*gas_gamma1;
}

template<typename REAL>
void
WxEulerEqn<REAL>::
RHS(unsigned N, REAL *geometry, REAL *normals, WxpDGGeometry<REAL> *quad, REAL *q, REAL *dq, REAL *rhs)
{

}

// instantiations
template class WxEulerEqn<float>;
template class WxEulerEqn<double>;
