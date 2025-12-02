// WarpX lib includes
#include <wxlogger.h>

// WarpX hyperbolic solver includes
#include "wxphmaxwelleqn.h"

// std includes
#include <cmath>

template<typename REAL>
void
WxPHMaxwellEqn<REAL>::
setup(const WxCryptSet& wxc)
{
  // light speed
  _c0 = wxc.template get<REAL>("c0");

  // error propagation speed for magnetic field
  if (wxc.has("gamma"))
    _gamma = wxc.template get<REAL>("gamma");
  else
    _gamma = 0.0;

  // error propagation speed for electric field
  if (wxc.has("chi"))
    _chi = wxc.template get<REAL>("chi");
  else
    _chi = 0.0;
}

template<typename REAL>
void
WxPHMaxwellEqn<REAL>::
rp(unsigned d, REAL *ql, REAL *qr, REAL *qauxl, REAL *qauxr, REAL *df, REAL **wave, REAL *s, REAL *amdq, REAL *apdq)
{
  unsigned iex=0, iey=0, iez=0, ibx=0, iby=0, ibz=0;
  REAL a1,a2,a3,a4,a5,a6,a7,a8;
  REAL delta[8];

  REAL c0 = _c0;
  REAL g = _gamma;
  REAL k = _chi;

  // set indexes into conserved variables array to correctly handle
  // X-Y-Z direction Reimann problems
  if (d==0) 
  {
    // X direction Reimann problem
    iex = 0;
    iey = 1;
    iez = 2;
    ibx = 3;
    iby = 4;
    ibz = 5;
  }
  else if (d==1)
  {
    // Y direction Reimann problem
    iex = 1;
    iey = 2;
    iez = 0;
    ibx = 4;
    iby = 5;
    ibz = 3;
  }
  else if (d==2)
  {
    // Z direction Reimann problem
    iex = 2;
    iey = 0;
    iez = 1;
    ibx = 5;
    iby = 3;
    ibz = 4;
  }

  // compute coefficients of the 4 eigenvectors
  delta[0] = df[iex];
  delta[1] = df[iey];
  delta[2] = df[iez];
  delta[3] = df[ibx];
  delta[4] = df[iby];
  delta[5] = df[ibz];
  delta[6] = df[6];
  delta[7] = df[7];

  a1 = 0.5*(delta[1]/c0+delta[5]);
  a2 = 0.5*(-delta[2]/c0+delta[4]);
  a3 = 0.5*(-delta[1]/c0+delta[5]);
  a4 = 0.5*(delta[2]/c0+delta[4]);
  a5 = 0.5*(delta[3]*c0+delta[7]);
  a6 = 0.5*(-delta[3]*c0+delta[7]);
  a7 = 0.5*(delta[0]/c0+delta[6]);
  a8 = 0.5*(-delta[0]/c0+delta[6]);

  // compute waves.
  //  there are six wave speeds and hence six waves.

  // Wave 1 corresponds to +c

  wave[iex][0] = 0;
  wave[iey][0] = c0*a1;
  wave[iez][0] = -c0*a2;
  wave[ibx][0] = 0;
  wave[iby][0] = a2;
  wave[ibz][0] = a1;
  wave[6][0]   = 0;
  wave[7][0]   = 0;
  s[0] = c0;

  // Wave 2 corresponds to -c

  wave[iex][1] = 0;
  wave[iey][1] = -c0*a3;
  wave[iez][1] = c0*a4;
  wave[ibx][1] = 0;
  wave[iby][1] = a4;
  wave[ibz][1] = a3;
  wave[6][1]   = 0;
  wave[7][1]   = 0;
  s[1] = -c0;

  // Wave 3 corresponds to c*g

  wave[iex][2] = 0;
  wave[iey][2] = 0;
  wave[iez][2] = 0;
  wave[ibx][2] = a5/c0;
  wave[iby][2] = 0;
  wave[ibz][2] = 0;
  wave[6][2]   = 0;
  wave[7][2]   = a5;
  s[2] = c0*g;

  // Wave 4 corresponds to -c*g

  wave[iex][3] = 0;
  wave[iey][3] = 0;
  wave[iez][3] = 0;
  wave[ibx][3] = -a6/c0;
  wave[iby][3] = 0;
  wave[ibz][3] = 0;
  wave[6][3] = 0;
  wave[7][3] = a6;
  s[3] = -c0*g;

  // Wave 5 corresponds to c*k

  wave[iex][4] = a7*c0;
  wave[iey][4] = 0;
  wave[iez][4] = 0;
  wave[ibx][4] = 0;
  wave[iby][4] = 0;
  wave[ibz][4] = 0;
  wave[6][4] = a7;
  wave[7][4] = 0;
  s[4] = c0*k;

  // Wave 6 corresponds to -c*k

  wave[iex][5] = -c0*a8;
  wave[iey][5] = 0;
  wave[iez][5] = 0;
  wave[ibx][5] = 0;
  wave[iby][5] = 0;
  wave[ibz][5] = 0;
  wave[6][5] = a8;
  wave[7][5] = 0;
  s[5] = -c0*k;

  // compute fluctuations
  for (unsigned m=0; m<8; ++m)
  {
    amdq[m] = 0.0; apdq[m] = 0.0;
    for (unsigned mw=0; mw<6; ++mw)
    {
      if (s[mw] < 0.0)
        // left going wave
        amdq[m] += s[mw]*wave[m][mw];
      else
        // right going wave
        apdq[m] += s[mw]*wave[m][mw];
    }
  }    
}

template<typename REAL>
void
WxPHMaxwellEqn<REAL>::
rpt(unsigned td, unsigned d, REAL *ql, REAL* qr,
  REAL *amdq, REAL* bmamdq, REAL* bpamdq,
  REAL *apdq, REAL* bmapdq, REAL* bpapdq)
{
  unsigned iex=0, iey=0, iez=0, ibx=0, iby=0, ibz=0;
  REAL a1,a2,a3,a4,a5,a6,a7,a8;
  REAL delta[8];

  REAL c0 = _c0;
  REAL g = _gamma;
  REAL k = _chi;

  // set indexes into conserved variables array to correctly handle
  // X-Y-Z direction Reimann problems
  if (d==0)
  { // x-direction Riemann problem
	  if (td==1)
	  {
		  iex = 1;
		  iey = 2;
		  iez = 0;
		  ibx = 4;
		  iby = 5;
		  ibz = 3;
	  }
	  else if (td==2)
	  {
		  iex = 2;
		  iey = 0;
		  iez = 1;
		  ibx = 5;
		  iby = 3;
		  ibz = 4;
	  }
	  else
			throw WxExcept("WxPhMaxwellEqn::rpt: transverse direction code not valid");

  }
  else if(d==1)
  { // y-direction Riemann problem
	  if (td==2)
	  {
		  iex = 2;
		  iey = 0;
		  iez = 1;
		  ibx = 5;
		  iby = 3;
		  ibz = 4;
	  }
	  else if (td==0)
	  {
		  iex = 0;
		  iey = 1;
		  iez = 2;
		  ibx = 3;
		  iby = 4;
		  ibz = 5;
	  }
	  else
			throw WxExcept("WxPhMaxwellEqn::rpt: transverse direction code not valid");
  }

  else if (d==2)
  { // z-direction Riemann problem
    if (td==0)
    {
  	  iex = 0;
  	  iey = 1;
  	  iez = 2;
  	  ibx = 3;
  	  iby = 4;
  	  ibz = 5;
    }
    else if (td==1)
    {
  	  iex = 1;
  	  iey = 2;
  	  iez = 0;
  	  ibx = 4;
  	  iby = 5;
  	  ibz = 3;
    }
    else
    	throw WxExcept("WxPhMaxwellEqn::rpt: transverse direction code not valid");
  }

  //
  // Split amdq
  //
  delta[0] = amdq[iex];
  delta[1] = amdq[iey];
  delta[2] = amdq[iez];
  delta[3] = amdq[ibx];
  delta[4] = amdq[iby];
  delta[5] = amdq[ibz];
  delta[6] = amdq[6];
  delta[7] = amdq[7];

  a1 = 0.5*(delta[1]/c0+delta[5]);
  a2 = 0.5*(-delta[2]/c0+delta[4]);
  a3 = 0.5*(-delta[1]/c0+delta[5]);
  a4 = 0.5*(delta[2]/c0+delta[4]);
  a5 = 0.5*(delta[3]*c0+delta[7]);
  a6 = 0.5*(-delta[3]*c0+delta[7]);
  a7 = 0.5*(delta[0]/c0+delta[6]);
  a8 = 0.5*(-delta[0]/c0+delta[6]);

  // compute up going fluctuation
  bpamdq[iex] = c0*k*c0*a7;
  bpamdq[iey] = c0*c0*a1;
  bpamdq[iez] = -c0*c0*a2;
  bpamdq[ibx] = g*a5;
  bpamdq[iby] = c0*a2;
  bpamdq[ibz] = c0*a1;
  bpamdq[6] = c0*k*a7;
  bpamdq[7] = c0*g*a5;

  // compute down going fluctuation
  bmamdq[iex] = c0*k*c0*a8;
  bmamdq[iey] = c0*c0*a3;
  bmamdq[iez] = -c0*c0*a4;
  bmamdq[ibx] = g*a6;
  bmamdq[iby] = -c0*a4;
  bmamdq[ibz] = -c0*a3;
  bmamdq[6] = -c0*k*a8;
  bmamdq[7] = -c0*g*a6;

  //
  // Split apdq
  //
  delta[0] = apdq[iex];
  delta[1] = apdq[iey];
  delta[2] = apdq[iez];
  delta[3] = apdq[ibx];
  delta[4] = apdq[iby];
  delta[5] = apdq[ibz];
  delta[6] = apdq[6];
  delta[7] = apdq[7];

  a1 = 0.5*(delta[1]/c0+delta[5]);
  a2 = 0.5*(-delta[2]/c0+delta[4]);
  a3 = 0.5*(-delta[1]/c0+delta[5]);
  a4 = 0.5*(delta[2]/c0+delta[4]);
  a5 = 0.5*(delta[3]*c0+delta[7]);
  a6 = 0.5*(-delta[3]*c0+delta[7]);
  a7 = 0.5*(delta[0]/c0+delta[6]);
  a8 = 0.5*(-delta[0]/c0+delta[6]);

  // compute up going fluctuation
  bpapdq[iex] = c0*k*c0*a7;
  bpapdq[iey] = c0*c0*a1;
  bpapdq[iez] = -c0*c0*a2;
  bpapdq[ibx] = g*a5;
  bpapdq[iby] = c0*a2;
  bpapdq[ibz] = c0*a1;
  bpapdq[6] = c0*k*a7;
  bpapdq[7] = c0*g*a5;

  // compute down going fluctuation
  bmapdq[iex] = c0*k*c0*a8;
  bmapdq[iey] = c0*c0*a3;
  bmapdq[iez] = -c0*c0*a4;
  bmapdq[ibx] = g*a6;
  bmapdq[iby] = -c0*a4;
  bmapdq[ibz] = -c0*a3;
  bmapdq[6] = -c0*k*a8;
  bmapdq[7] = -c0*g*a6;
}

template<typename REAL>
void
WxPHMaxwellEqn<REAL>::
rptc(unsigned td, unsigned d, REAL *ql, REAL* qr,REAL *soc, REAL* bms, REAL* bps)
{
	  unsigned iex=0, iey=0, iez=0, ibx=0, iby=0, ibz=0;
	  REAL a1,a2,a3,a4,a5,a6,a7,a8;
	  REAL delta[8];

	  REAL c0 = _c0;
	  REAL g = _gamma;
	  REAL k = _chi;

	  // set indexes into conserved variables array to correctly handle
	  // X-Y-Z direction Reimann problems
	  if (td==0)
	  {
		  iex = 0;
		  iey = 1;
		  iez = 2;
		  ibx = 3;
		  iby = 4;
		  ibz = 5;
	  }
	  else if (td==1)
	  {
		  iex = 1;
		  iey = 2;
		  iez = 0;
		  ibx = 4;
		  iby = 5;
		  ibz = 3;
	  }
	  else if (td==2)
	  {
		  iex = 2;
		  iey = 0;
		  iez = 1;
		  ibx = 5;
		  iby = 3;
		  ibz = 4;
	  }
	  else
		throw WxExcept("WxEulerEqn::rptc: transverse direction code not valid");

	  //
	  // Split amdq
	  //
	  delta[0] = soc[iex];
	  delta[1] = soc[iey];
	  delta[2] = soc[iez];
	  delta[3] = soc[ibx];
	  delta[4] = soc[iby];
	  delta[5] = soc[ibz];
	  delta[6] = soc[6];
	  delta[7] = soc[7];

	  a4 = 0.5*(delta[2]/c0+delta[4]);
	  a1 = 0.5*(delta[1]/c0+delta[5]);
	  a5 = 0.5*(delta[3]*c0+delta[7]);
	  a7 = 0.5*(delta[0]/c0+delta[6]);

	  a2 = 0.5*(-delta[2]/c0+delta[4]);
	  a3 = 0.5*(-delta[1]/c0+delta[5]);
	  a6 = 0.5*(-delta[3]*c0+delta[7]);
	  a8 = 0.5*(-delta[0]/c0+delta[6]);

	  // compute up going fluctuation
	  bps[iex] = c0*k*c0*a7;
	  bps[iey] = c0*c0*a1;
	  bps[iez] = -c0*c0*a2;
	  bps[ibx] = g*a5;
	  bps[iby] = c0*a2;
	  bps[ibz] = c0*a1;
	  bps[6] = c0*k*a7;
	  bps[7] = c0*g*a5;

	  // compute down going fluctuation
	  bms[iex] = c0*k*c0*a8;
	  bms[iey] = c0*c0*a3;
	  bms[iez] = -c0*c0*a4;
	  bms[ibx] = g*a6;
	  bms[iby] = -c0*a4;
	  bms[ibz] = -c0*a3;
	  bms[6] = -c0*k*a8;
	  bms[7] = -c0*g*a6;
}

template<typename REAL>
void
WxPHMaxwellEqn<REAL>::
flux(unsigned d, REAL *x, REAL *q, REAL *qaux, REAL *f)
{
  unsigned iex=0, iey=0, iez=0, ibx=0, iby=0, ibz=0;

  // set indexes into conserved variables array to correctly handle
  // X-Y-Z direction Reimann problems
  if (d==0) 
  {
    // X direction flux
    iex = 0;
    iey = 1;
    iez = 2;
    ibx = 3;
    iby = 4;
    ibz = 5;
  }
  else if (d==1)
  {
    // Y direction flux
    iex = 1;
    iey = 2;
    iez = 0;
    ibx = 4;
    iby = 5;
    ibz = 3;
  }
  else if (d==2)
  {
    // Z direction flux
    iex = 2;
    iey = 0;
    iez = 1;
    ibx = 5;
    iby = 3;
    ibz = 4;
  }

//  REAL AA = q[ibx];
//  REAL AB = q[7];

  f[iex] = _chi*_c0*_c0*q[6];
  f[iey] = _c0*_c0*q[ibz];
  f[iez] = -_c0*_c0*q[iby];
  f[ibx] = _gamma*q[7];
  f[iby] = -q[iez];
  f[ibz] = q[iey];
  f[6] = _chi*q[iex];
  f[7] = _gamma*_c0*_c0*q[ibx];
}

template<typename REAL>
void
WxPHMaxwellEqn<REAL>::
DGnumericalFlux(REAL *normals, REAL *qM, REAL *qP, REAL *nflux, REAL *maxSpeed)
{
    REAL *xc, *qaux;
    REAL fM[8], fP[8], gM[8], gP[8]; // x/y Fluxes

    // evaluate fluxes
    this->flux(0, xc, qM, qaux, fM);
    this->flux(0, xc, qP, qaux, fP);
    this->flux(1, xc, qM, qaux, gM);
    this->flux(1, xc, qP, qaux, gP);

    // determine the fastest propagating wave speed
    REAL lambda = dmax(_chi*_c0,_gamma*_c0,_c0);

    // Lax-Frederick fluxes
    for(unsigned comp=0; comp<meqn(); comp++)
        nflux[comp] = 0.5*(normals[0]*(fM[comp]+fP[comp]) + normals[1]*(gM[comp]+gP[comp]) + lambda*(qM[comp]-qP[comp]));

    *maxSpeed = lambda;
}

template<typename REAL>
void
WxPHMaxwellEqn<REAL>::
eigenSystem(unsigned d, REAL *q, REAL *ev, REAL **lev, REAL **rev)
{
  unsigned iex=0, iey=0, iez=0, ibx=0, iby=0, ibz=0;
  REAL gamma, kappa;

  // set indexes into conserved variables array to correctly handle
  // X-Y-Z direction Reimann problems
  if (d==0) 
  {
    // X direction Reimann problem
    iex = 0;
    iey = 1;
    iez = 2;
    ibx = 3;
    iby = 4;
    ibz = 5;
  }
  else if (d==1)
  {
    // Y direction Reimann problem
    iex = 1;
    iey = 2;
    iez = 0;
    ibx = 4;
    iby = 5;
    ibz = 3;
  }
  else if (d==2)
  {
    // Z direction Reimann problem
    iex = 2;
    iey = 0;
    iez = 1;
    ibx = 5;
    iby = 3;
    ibz = 4;
  }

  gamma = q[6];
  kappa = q[7];

  // eigenvalues
  ev[0] = ev[1] = _c0;
  ev[2] = ev[3] = -_c0;
  ev[4] = _c0*gamma;
  ev[5] = -_c0*gamma;
  ev[6] = _c0*kappa;
  ev[7] = -_c0*kappa;

  // right eigenvectors: these are arranged as column vectors
  rev[iex][0]  = 0.;
  rev[iey][0] = _c0;
  rev[iez][0] = 0.;
  rev[ibx][0] = 0.;
  rev[iby][0] = 0.;
  rev[ibz][0] = 1.;
  rev[6][0] = 0.;
  rev[7][0] = 0.;

  rev[iex][1]  = 0.;
  rev[iey][1] = 0.;
  rev[iez][1] = -_c0;
  rev[ibx][1] = 0.;
  rev[iby][1] = 1.;
  rev[ibz][1] = 0.;
  rev[6][1] = 0.;
  rev[7][1] = 0.;

  rev[iex][2]  = 0.;
  rev[iey][2] = -_c0;
  rev[iez][2] = 0.;
  rev[ibx][2] = 0.;
  rev[iby][2] = 0.;
  rev[ibz][2] = 1.;
  rev[6][2] = 0.;
  rev[7][2] = 0.;

  rev[iex][3]  = 0.;
  rev[iey][3] = 0.;
  rev[iez][3] = _c0;
  rev[ibx][3] = 0.;
  rev[iby][3] = 1.;
  rev[ibz][3] = 0.;
  rev[6][3] = 0.;
  rev[7][3] = 0.;

  rev[iex][4]  = 0.;
  rev[iey][4] = 0.;
  rev[iez][4] = 0.;
  rev[ibx][4] = 1/_c0;
  rev[iby][4] = 0.;
  rev[ibz][4] = 0.;
  rev[6][4] = 0.;
  rev[7][4] = 1.;

  rev[iex][5]  = 0.;
  rev[iey][5] = 0.;
  rev[iez][5] = 0.;
  rev[ibx][5] = -1/_c0;
  rev[iby][5] = 0.;
  rev[ibz][5] = 0.;
  rev[6][5] = 0.;
  rev[7][5] = 1.;

  rev[iex][6]  = _c0;
  rev[iey][6] = 0.;
  rev[iez][6] = 0.;
  rev[ibx][6] = 0.;
  rev[iby][6] = 0.;
  rev[ibz][6] = 0.;
  rev[6][6] = 1.;
  rev[7][6] = 0.;

  rev[iex][7]  = -_c0;
  rev[iey][7] = 0.;
  rev[iez][7] = 0.;
  rev[ibx][7] = 0.;
  rev[iby][7] = 0.;
  rev[ibz][7] = 0.;
  rev[6][7] = 1.;
  rev[7][7] = 0.;

  // left eigenvectors: these are arranged as row vectors
  lev[0][iex]  = 0.;
  lev[0][iey] = 1/(2.*_c0);
  lev[0][iez] = 0.;
  lev[0][ibx] = 0.;
  lev[0][iby] = 0.;
  lev[0][ibz] = 1./2.;
  lev[0][6] = 0.;
  lev[0][7] = 0.;

  lev[1][iex]  = 0.;
  lev[1][iey] = 0.;
  lev[1][iez] = -1/(2.*_c0);
  lev[1][ibx] = 0.;
  lev[1][iby] = 1./2.;
  lev[1][ibz] = 0.;
  lev[1][6] = 0.;
  lev[1][7] = 0.;

  lev[2][iex]  = 0.;
  lev[2][iey] = -1/(2.*_c0);
  lev[2][iez] = 0.;
  lev[2][ibx] = 0.;
  lev[2][iby] = 0.;
  lev[2][ibz] = 1./2.;
  lev[2][6] = 0.;
  lev[2][7] = 0.;

  lev[3][iex]  = 0.;
  lev[3][iey] = 0.;
  lev[3][iez] = 1/(2.*_c0);
  lev[3][ibx] = 0.;
  lev[3][iby] = 1./2.;
  lev[3][ibz] = 0.;
  lev[3][6] = 0.;
  lev[3][7] = 0.;

  lev[4][iex]  = 0.;
  lev[4][iey] = 0.;
  lev[4][iez] = 0.;
  lev[4][ibx] = _c0/2.;
  lev[4][iby] = 0.;
  lev[4][ibz] = 0.;
  lev[4][6] = 0.;
  lev[4][7] = 1./2.;

  lev[5][iex]  = 0.;
  lev[5][iey] = 0.;
  lev[5][iez] = 0.;
  lev[5][ibx] = -_c0/2.;
  lev[5][iby] = 0.;
  lev[5][ibz] = 0.;
  lev[5][6] = 0.;
  lev[5][7] = 1./2.;

  lev[6][iex]  = 1/(2.*_c0);
  lev[6][iey] = 0.;
  lev[6][iez] = 0.;
  lev[6][ibx] = 0.;
  lev[6][iby] = 0.;
  lev[6][ibz] = 0.;
  lev[6][6] = 1./2.;
  lev[6][7] = 0.;

  lev[7][iex]  = -1/(2.*_c0);
  lev[7][iey] = 0.;
  lev[7][iez] = 0.;
  lev[7][ibx] = 0.;
  lev[7][iby] = 0.;
  lev[7][ibz] = 0.;
  lev[7][6] = 1./2.;
  lev[7][7] = 0.;

}

template<typename REAL>
void
WxPHMaxwellEqn<REAL>::
primitiveVariables(REAL *qCons, REAL *qPrim)
{
    qPrim[0] = qCons[0];
    qPrim[1] = qCons[1];
    qPrim[2] = qCons[2];
    qPrim[3] = qCons[3];
    qPrim[4] = qCons[4];
    qPrim[5] = qCons[5];
    qPrim[6] = qCons[6];
    qPrim[7] = qCons[7];
}

template<typename REAL>
void
WxPHMaxwellEqn<REAL>::
limiterTuAndAliabadi(REAL *avgCons, REAL *avgPrim, REAL *dGrads, REAL *limitedValues)
{

    limitedValues[0] = avgPrim[0];
    limitedValues[1] = avgPrim[1];
    limitedValues[2] = avgPrim[2];
    limitedValues[3] = avgPrim[3];
    limitedValues[4] = avgPrim[4];
    limitedValues[5] = avgPrim[5];
    limitedValues[6] = avgPrim[6];
    limitedValues[7] = avgPrim[7];

}

// instantiations
//template class WxPHMaxwellEqn<float>;
template class WxPHMaxwellEqn<double>;
