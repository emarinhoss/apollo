// WarpX lib includes
#include <wxlogger.h>
#include <wxmath.h>

// WarpX hyperbolic solver includes
#include "wxadvectionunstuctured.h"

// std includes
#include <cmath>

template<typename REAL>
void
WxAdvectionUnstructuredEqn<REAL>::
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

template <typename REAL>
void
WxAdvectionUnstructuredEqn<REAL>::rotateToLocalFrame(REAL norm[3], REAL tan1[3],
  REAL tan2[3], REAL *vin, REAL *vout)
{
//	REAL _uout[3];
//
//	_uout[0] = _ux*norm[0] + _uy*norm[1] + _uz*norm[2];
//	_uout[1] = _ux*tan1[0] + _uy*tan1[1] + _uz*tan1[2];
//	_uout[2] = _ux*tan2[0] + _uy*tan2[1] + _uz*tan2[2];

//	if (abs(norm[1])==1)
//		std::cout << "u_vec_in = (" << _ux << "," << _uy << "," << _uz <<")" << std::endl;
//
//	_ux = _uout[0];
//	_uy = _uout[1];
//	_uz = _uout[2];

//	if (abs(norm[1])==1)
//		std::cout << "u_vec_out = (" << _ux << "," << _uy << "," << _uz <<")" << std::endl;

    vout[0] = vin[0];

}

template <typename REAL>
void
WxAdvectionUnstructuredEqn<REAL>::rotateToGlobalFrame(REAL norm[3], REAL tan1[3],
  REAL tan2[3], REAL *vin, REAL *vout)
{
//	  REAL t1x = tan1[0];
//	  REAL t1y = tan1[1];
//	  REAL t1z = tan1[2];
//
//	  REAL t2x = tan2[0];
//	  REAL t2y = tan2[1];
//	  REAL t2z = tan2[2];
//
//	  REAL nx = norm[0];
//	  REAL ny = norm[1];
//	  REAL nz = norm[2];
//
//	  REAL   Num11 = t1y*t2z-t2y*t1z;
//	  REAL   Num12 = t2x*t1z-t1x*t2z;
//	  REAL   Num13 = t1x*t2y-t1y*t2x;
//	  REAL   Num21 = nz*t2y-ny*t2z;
//	  REAL   Num22 = nx*t2z-nz*t2x;
//	  REAL   Num23 = ny*t2x-nx*t2y;
//	  REAL   Num31 = ny*t1z-nz*t1y;
//	  REAL   Num32 = nz*t1x-nx*t1z;
//	  REAL   Num33 = nx*t1y-ny*t1x;
//
//	  REAL   Dem = nx*Num11+ny*Num12+nz*Num13;
//
//	  REAL _uxout = (_ux*Num11 + _uy*Num12 + _uz*Num13)/Dem;
//	  REAL _uyout = (_ux*Num21 + _uy*Num22 + _uz*Num23)/Dem;
//	  REAL _uzout = (_ux*Num31 + _uy*Num32 + _uz*Num33)/Dem;
//
//	  _ux = _uxout;
//	  _uy = _uyout;
//	  _uz = _uzout;

      vout[0] = vin[0];
}

template <typename REAL>
void
WxAdvectionUnstructuredEqn<REAL>::rotateToLocalFrame(REAL norm[3], REAL *vin, REAL *vout)
{
  // rotation function to rotate 'vin' to local coordinate system
  vout[0] = vin[0]*norm[0];// + vin[0]*norm[1];

}

template <typename REAL>
void
WxAdvectionUnstructuredEqn<REAL>::rotateToGlobalFrame(REAL norm[3], REAL *vin, REAL *vout)
{
  vout[0] = vin[0]*norm[0];// - vin[0]*norm[1];
}

template<typename REAL>
void
WxAdvectionUnstructuredEqn<REAL>::
rp(unsigned d, REAL *ql, REAL *qr,
        REAL *qauxl, REAL *qauxr, REAL *df,
        REAL **wave, REAL *s, REAL *amdq, REAL *apdq)
{
  unsigned meqn = this->meqn();
  unsigned mwave = this->mwave();
  REAL norm[3]; s[0]=0.0;

  for(unsigned kk=0; kk<3; kk++)
      norm[kk] = qauxl[kk];
  // compute waves

  // wave 1
  wave[0][0] = df[0];

  s[0] = _ux*norm[0] + _uy*norm[1] + _uz*norm[2];

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
WxAdvectionUnstructuredEqn<REAL>::
rpt(unsigned td, unsigned d, REAL *ql, REAL* qr,
  REAL *amdq, REAL* bmamdq, REAL* bpamdq,
  REAL *apdq, REAL* bmapdq, REAL* bpapdq)
{
    unsigned meqn = this->meqn();
//	unsigned mwave = this->mwave();

    REAL norm[3],tan1[3],tan2[3];

    this->getFaceVectors(norm, tan1, tan2);

    REAL uout[3], vtrans;

    //bool xface = ((fabs(norm[0])==1.0)&&((norm[1])==0.0)&&((norm[2])==0.0));
    //bool yface = ((fabs(norm[1])==1.0)&&((norm[0])==0.0)&&((norm[2])==0.0));
    //bool zface = ((fabs(norm[2])==1.0)&&((norm[1])==0.0)&&((norm[0])==0.0));
    //bool faceisnotcartesian = ((xface!=1)&&(yface!=1)&&(zface!=1));
    //bool havefluctuation = (apdq[0]!=0.0)||(amdq[0]!=0.0);

//	if (d==0)
//		std::cout << "d = " << d << ",\t norm = (" << norm[0] << ", " << norm[1] << "," << norm[2] << "), tan1 = (" << tan1[0] << "," << tan1[1] << "," << tan1[2] << "), tan2 = (" << tan2[0] << "," << tan2[1] << ","<<  tan2[2] << ")" << std::endl;
//
//
//	if (havefluctuation&&faceisnotcartesian)
//		std::cout << "d = " << d << ",\t norm = (" << norm[0] << ", " << norm[1] << "," << norm[2] << "), tan1 = (" << tan1[0] << "," << tan1[1] << "," << tan1[2] << "), tan2 = (" << tan2[0] << "," << tan2[1] << ","<<  tan2[2] << ")" << std::endl;

    uout[0] = _ux*norm[0] + _uy*norm[1] + _uz*norm[2];
    uout[1] = _ux*tan1[0] + _uy*tan1[1] + _uz*tan1[2];
    uout[2] = _ux*tan2[0] + _uy*tan2[1] + _uz*tan2[2];

//    if (d==0)
//    {
//    	if (td==1)
//    		vtrans = uout[1];
//    	else if (td==2)
//    		vtrans = uout[2];
//    	else
//    	{
//    		std::cout << "(d,td) = (" << d << "," << td << ")" << std::endl;
//    		throw WxExcept("WxAdvectionGenGeomEqn::rpt : invalid transverse direction\n");
//    	}
//    }
//    else if (d==1)
//    {
//    	if (td==2)
//    		vtrans = uout[1];
//    	else if (td==0)
//    		vtrans = uout[2];
//    	else
//    	{
//    		std::cout << "(d,td) = (" << d << "," << td << ")" << std::endl;
//    		throw WxExcept("WxAdvectionGenGeomEqn::rpt : invalid transverse direction\n");
//    	}
//    }
//    else if (d==2)
//    {
//    	if (td==0)
//    		vtrans = uout[1];
//    	else if (td==1)
//    		vtrans = uout[2];
//    	else
//    	{
//    		std::cout << "(d,td) = (" << d << "," << td << ")" << std::endl;
//    		throw WxExcept("WxAdvectionGenGeomEqn::rpt : invalid transverse direction\n");
//    	}
//    }

//    if (((apdq[0]!=0.0)||(amdq[0]!=0.0))&&(vtrans!=0))
//    	std::cout << "apdq = " << apdq[0]  << ", \t amdq = " << amdq[0] <<", \t u = (" << _ux << ", " << _uy << ", " << _uz << ") and u_local = ("<< uout[0] << ", "<< uout[1] << ", " << uout[2] << "), and finally s[0] = "<< vtrans << std::endl;

  vtrans = uout[0];

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
WxAdvectionUnstructuredEqn<REAL>::
rptc(unsigned td, unsigned d, REAL *ql, REAL* qr,
  REAL *s, REAL* bms, REAL* bps)
{
    REAL norm[3],tan1[3],tan2[3];

    this->getFaceVectors(norm, tan1, tan2);

    REAL uout[3], vtrans;

    uout[0] = _ux*norm[0] + _uy*norm[1] + _uz*norm[2];
    uout[1] = _ux*tan1[0] + _uy*tan1[1] + _uz*tan1[2];
    uout[2] = _ux*tan2[0] + _uy*tan2[1] + _uz*tan2[2];

//    if (d==0)
//    {
//    	if (td==1)
//    		vtrans = uout[1];
//    	else if (td==2)
//    		vtrans = uout[2];
//    	else
//    	{
//    		std::cout << "(d,td) = (" << d << "," << td << ")" << std::endl;
//    		throw WxExcept("WxAdvectionGenGeomEqn::rptc : invalid transverse direction\n");
//    	}
//    }
//    else if (d==1)
//    {
//    	if (td==2)
//    		vtrans = uout[1];
//    	else if (td==0)
//    		vtrans = uout[2];
//    	else
//    	{
//    		std::cout << "(d,td) = (" << d << "," << td << ")" << std::endl;
//    		throw WxExcept("WxAdvectionGenGeomEqn::rpt : invalid transverse direction\n");
//    	}
//    }
//    else if (d==2)
//    {
//    	if (td==0)
//    		vtrans = uout[1];
//    	else if (td==1)
//    		vtrans = uout[2];
//    	else
//    	{
//    		std::cout << "(d,td) = (" << d << "," << td << ")" << std::endl;
//    		throw WxExcept("WxAdvectionGenGeomEqn::rpt : invalid transverse direction\n");
//    	}
//    }

  vtrans = uout[0];

  REAL vtransm, vtransp;
  vtransm = dmin<REAL>(vtrans, 0.0);
  vtransp = dmax<REAL>(vtrans, 0.0);

  // split s
  bms[0] = vtransm*s[0];
  bps[0] = vtransp*s[0];

//  if ((bps[0]!=0.0)||(bms[0]!=0.0))
//	  std::cout << "vtrans = " << vtrans << ", \t s[0] = " << s[0] << ", \t bps[0] = " << bps[0] << ", \t bms[0] = " << bms[0] << std::endl;
}

template<typename REAL>
void
WxAdvectionUnstructuredEqn<REAL>::
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
WxAdvectionUnstructuredEqn<REAL>::
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
template class WxAdvectionUnstructuredEqn<float>;
template class WxAdvectionUnstructuredEqn<double>;
