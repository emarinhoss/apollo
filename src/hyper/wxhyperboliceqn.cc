// hyper includes
#include "wxhyperboliceqn.h"

template <typename REAL>
WxHyperbolicEqn<REAL>::WxHyperbolicEqn(const std::string& name)
  : _name(name) 
{
  for (unsigned i=0; i<3; ++i)
  {
    _xl[i] = _xr[i] = 0.0;
  }
}

template <typename REAL>
WxHyperbolicEqn<REAL>::~WxHyperbolicEqn()
{
}

template <typename REAL>
std::string
WxHyperbolicEqn<REAL>::name() const
{
  return _name;
}

template <typename REAL>
void
WxHyperbolicEqn<REAL>::rotateToLocalFrame(REAL norm[3], REAL tan1[3], REAL tan2[3], REAL *vin, REAL *vout)
{
  WxExcept wxe("WxWxHyperbolicEqn::rotateToLocalFrame: ");
  wxe << " Not implemented for equation system " << _name;
  throw wxe;
}

template <typename REAL>
void
WxHyperbolicEqn<REAL>::rotateToGlobalFrame(REAL norm[3], REAL tan1[3], REAL tan2[3], REAL *vin, REAL *vout)
{
  WxExcept wxe("WxWxHyperbolicEqn::rotateToGlobalFrame: ");
  wxe << " Not implemented for equation system " << _name;
  throw wxe;
}

template <typename REAL>
void
WxHyperbolicEqn<REAL>::rotateToLocalFrame(REAL norm[3], REAL *vin, REAL *vout)
{
  WxExcept wxe("WxWxHyperbolicEqn::rotateToLocalFrame: ");
  wxe << " Not implemented for equation system " << _name;
  throw wxe;
}

template <typename REAL>
void
WxHyperbolicEqn<REAL>::rotateToGlobalFrame(REAL norm[3], REAL *vin, REAL *vout)
{
  WxExcept wxe("WxWxHyperbolicEqn::rotateToGlobalFrame: ");
  wxe << " Not implemented for equation system " << _name;
  throw wxe;
}

template <typename REAL>
void
WxHyperbolicEqn<REAL>::setCoordLeftCell(unsigned ndim, REAL *x) 
{
  for (unsigned i=0; i<ndim; ++i)
    _xl[i] = x[i];
}

template <typename REAL>
void
WxHyperbolicEqn<REAL>::setCoordRightCell(unsigned ndim, REAL *x) 
{
  for (unsigned i=0; i<ndim; ++i)
    _xr[i] = x[i];
}

template <typename REAL>
void
WxHyperbolicEqn<REAL>::coordLeftCell(REAL *x)
{
  for (unsigned i=0; i<3; ++i)
    x[i] = _xl[i];
}

template <typename REAL>
void
WxHyperbolicEqn<REAL>::coordRightCell(REAL *x) 
{
  for (unsigned i=0; i<3; ++i)
    x[i] = _xr[i];
}

template <typename REAL>
void
WxHyperbolicEqn<REAL>::flux(unsigned d, REAL *x, REAL *q, REAL *qaux, REAL *f) 
{
  WxExcept wxe;
  wxe << "Flux function flux() of equation system" << _name
      << " not implemented." << std::endl;
  throw wxe;        
}

template <typename REAL>
void
WxHyperbolicEqn<REAL>::fluxJacobian(unsigned d, REAL *x, REAL *q, REAL *qaux, REAL **f)
{
  WxExcept wxe;
  wxe << "Flux Jacobian function fluxJacobian() of equation system" << _name
      << " not implemented." << std::endl;
  throw wxe;
}

template <typename REAL>
void
WxHyperbolicEqn<REAL>::edgefluxgengeom(unsigned d, REAL *x, REAL *q, REAL *qaux, REAL *f) 
{
  WxExcept wxe;
  wxe << "Edge flux gen geom function flux() of equation system" << _name
      << " not implemented." << std::endl;
  throw wxe;        
}

template <typename REAL>
void
WxHyperbolicEqn<REAL>::RHS(unsigned N, REAL *geometry, REAL *normals, WxpDGGeometry<REAL> *quad, REAL *q, REAL *dq, REAL *rhs)
{
  WxExcept wxe;
  wxe << "RHS calculations for equation system" << _name
      << " not implemented." << std::endl;
  throw wxe;
}

template <typename REAL>
void
WxHyperbolicEqn<REAL>::rpt(unsigned td, unsigned d, REAL *ql, REAL* qr,
      REAL *amdq, REAL* bmamdq, REAL* bpamdq,
      REAL *apdq, REAL* bmapdq, REAL* bpapdq) 
{
  // throw an exception if called and derived class does not define it
  WxExcept wxe;
  wxe << "Transverse solver function rpt() of equation system " << _name
      << " not implemented." << std::endl;
  throw wxe;        
}

template <typename REAL>
void
WxHyperbolicEqn<REAL>::rptc(unsigned td, unsigned d, REAL *ql, REAL* qr,
      REAL *s, REAL* bms, REAL* bps)
{
  // throw an exception if called and derived class does not define it
  WxExcept wxe;
  wxe << "Transverse correction solver function rptc() of equation system " << _name
      << " not implemented." << std::endl;
  throw wxe;
}

template <typename REAL>
void
WxHyperbolicEqn<REAL>::rplimit(unsigned d, REAL *ql, REAL* qr, REAL *ql1, REAL* qr1, 
                               REAL *qauxl, REAL *qauxr, REAL *df, 
                               REAL **wave, REAL *s, REAL *amdq, REAL *apdq)
{
  // throw an exception if called and derived class does not define it
  WxExcept wxe;
  wxe << "Limiter-applied rplimit() of equation system " << _name
      << " not implemented." << std::endl;
  throw wxe;        
}

template <typename REAL>
void
WxHyperbolicEqn<REAL>::eigenSystem(unsigned d, REAL *q, REAL *ev, REAL **lev, REAL **rev) 
{
  // throw an exception if called and derived class does not define it
  WxExcept wxe;
  wxe << "Eigensystem function eigenSystem() of equation system " << _name
      << " not implemented." << std::endl;
  throw wxe;
}

template <typename REAL>
void
WxHyperbolicEqn<REAL>::
setFaceVectors(REAL norm[3], REAL tan1[3], REAL tan2[3])
{

	_norm[0] = norm[0];
	_norm[1] = norm[1];
	_norm[2] = norm[2];

	_tan1[0] = tan1[0];
	_tan1[1] = tan1[1];
	_tan1[2] = tan1[2];

	_tan2[0] = tan2[0];
	_tan2[1] = tan2[1];
	_tan2[2] = tan2[2];
}

template <typename REAL>
void
WxHyperbolicEqn<REAL>::
getFaceVectors(REAL norm[3], REAL tan1[3], REAL tan2[3])
{

	norm[0] = _norm[0];
	norm[1] = _norm[1];
	norm[2] = _norm[2];

	tan1[0] = _tan1[0];
	tan1[1] = _tan1[1];
	tan1[2] = _tan1[2];

	tan2[0] = _tan2[0];
	tan2[1] = _tan2[1];
	tan2[2] = _tan2[2];
}

template <typename REAL>
void
WxHyperbolicEqn<REAL>::
setIndices(int idx[3])
{
	_idx[0] = idx[0];
	_idx[1] = idx[1];
	_idx[2] = idx[2];
}

template <typename REAL>
void
WxHyperbolicEqn<REAL>::
getIndices(int idx[3])
{
	idx[0] = _idx[0];
	idx[1] = _idx[1];
	idx[2] = _idx[2];
}
// instantiations
template class WxHyperbolicEqn<float>;
template class WxHyperbolicEqn<double>;
