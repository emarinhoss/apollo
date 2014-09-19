// WarpX lib includes
#include <wxlogger.h>
#include <wxmath.h>

// WarpX hyperbolic solver includes
#include "wxadvectioneqn.h"

// std includes
#include <cmath>

template<typename REAL>
void
WxAdvectionEqn<REAL>::
setup(const WxCryptSet& wxc)
{
  if (wxc.has("ux"))
    _ux = wxc.template get<REAL>("ux");
  else
    _ux = 0.0;

  if (wxc.has("uy"))
    _uy = wxc.template get<REAL>("uy");
  else
    _uy = 0.0;

  if (wxc.has("uz"))
    _uz = wxc.template get<REAL>("uz");
  else
    _uz = 0.0;
}

template<typename REAL>
void
WxAdvectionEqn<REAL>::
rp(unsigned d, REAL *ql, REAL *qr, 
		REAL *qauxl, REAL *qauxr, REAL *df, 
		REAL **wave, REAL *s, REAL *amdq, REAL *apdq)
{
  unsigned meqn = this->meqn();
  unsigned mwave = this->mwave();
  // compute waves

  // wave 1 
  wave[0][0] = df[0];
  if (d == 0)
    s[0] = _ux;
  else if (d == 1)
    s[0] = _uy;
  else
    s[0] = _uz;

  // compute fluctuations
  for (unsigned m=0; m<meqn; ++m)
  {
    amdq[m] = 0.0; apdq[m] = 0.0;
    for (unsigned mw=0; mw<mwave; ++mw)
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
WxAdvectionEqn<REAL>::
rplimit(unsigned d, REAL *ql, REAL *qr, REAL *ql1, REAL *qr1, REAL *qauxl, REAL *qauxr, REAL *df, REAL **wave, REAL *s, REAL *amdq, REAL *apdq)
{
  unsigned meqn = this->meqn();
  unsigned mwave = this->mwave();
  // compute waves

  // wave 1 
  wave[0][0] = df[0];
  if (d == 0)
    s[0] = _ux;
  else if (d == 1)
    s[0] = _uy;
  else
    s[0] = _uz;

  // compute fluctuations
  for (unsigned m=0; m<meqn; ++m)
  {
    amdq[m] = 0.0; apdq[m] = 0.0;
    for (unsigned mw=0; mw<mwave; ++mw)
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
WxAdvectionEqn<REAL>::
rpt(unsigned td, unsigned d, REAL *ql, REAL* qr,
  REAL *amdq, REAL* bmamdq, REAL* bpamdq,
  REAL *apdq, REAL* bmapdq, REAL* bpapdq)
{
	unsigned meqn = this->meqn();
//	unsigned mwave = this->mwave();


  REAL vtrans = _ux;
  if (td == 0)
    // use X-direction speed
    vtrans = _ux;
  else if (td == 1)
    // use Y-direction speed
    vtrans = _uy;
  else
    // use Z-direction speed
    vtrans = _uz;

  for (unsigned m=0; m<meqn; ++m)
  {
    bpapdq[m] = 0.0;
    bmapdq[m] = 0.0;
    bpamdq[m] = 0.0;
    bmamdq[m] = 0.0;

    if (vtrans < 0.0)
    {
    //  left going wave
	    bmamdq[m] += vtrans*amdq[m];
	    bmapdq[m] += vtrans*apdq[m];
    }
    else
    {
    //  right going wave
	    bpamdq[m] += vtrans*amdq[m];
	    bpapdq[m] += vtrans*apdq[m];
    }
  }
}


template<typename REAL>
void
WxAdvectionEqn<REAL>::
rptc(unsigned td, unsigned d, REAL *ql, REAL* qr,
  REAL *s, REAL* bms, REAL* bps)
{
  REAL vtrans = _ux;
  if (d == 0)
    // use X-direction speed
    vtrans = _ux;
  else if (d == 1)
    // use Y-direction speed
    vtrans = _uy;
  else
    // use Z-direction speed
    vtrans = _uz;

  REAL vtransm, vtransp;
  vtransm = dmin<REAL>(vtrans, 0.0);
  vtransp = dmax<REAL>(vtrans, 0.0);

  // split s
  bms[0] = vtransm*s[0];
  bps[0] = vtransp*s[0];
}

template<typename REAL>
void
WxAdvectionEqn<REAL>::
flux(unsigned d, REAL *x, REAL *q, REAL *qaux, REAL *f)
{
  if (d==0)
    f[0] = _ux*q[0];
  else if (d==1)
    f[0] = _uy*q[0];
  else
    f[0] = _uz*q[0];
}

template<typename REAL>
void
WxAdvectionEqn<REAL>::
fluxJacobian(unsigned d, REAL *x, REAL *q, REAL *qaux, REAL **f)
{
    if (d==0)
      f[0][0] = _ux;
    else if (d==1)
      f[0][0] = _uy;
    else
      f[0][0] = _uz;
}

template<typename REAL>
void
WxAdvectionEqn<REAL>::
eigenSystem(unsigned d, REAL *q, REAL *ev, REAL **lev, REAL **rev)
{
  // eigenvalues
  if (d==0)
    ev[0] = _ux;
  else if (d==1)
    ev[0] = _uy;
  else 
    ev[0] = _uz;

  // right eigenvectors: these are arranged as column vectors
  rev[0][0]  = 1.;

  // left eigenvectors: these are arranged as row vectors
  lev[0][0]  = 1.;
} 

// instantiations
template class WxAdvectionEqn<float>;
template class WxAdvectionEqn<double>;
