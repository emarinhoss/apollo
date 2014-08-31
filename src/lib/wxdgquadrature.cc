#include "wxdgquadrature.h"

// WarpX lib includes
#include <wxmath.h>

template <typename REAL>
WxDGQuadrature<REAL>::WxDGQuadrature(unsigned meqn, unsigned spatialOrder, bool useGaussian)
  : _meqn(meqn), _spatialOrder(spatialOrder),
    _ind1D(WxRange(0,meqn, 0,(_spatialOrder+1))),
    _ind2D(WxRange(0,meqn, 0,(_spatialOrder+1), 0,(_spatialOrder+1))),
    _ind3D(WxRange(0,meqn, 0,(_spatialOrder+1), 0,(_spatialOrder+1), 0,(_spatialOrder+1)))
{
  _numCoeffPerElem = _spatialOrder + 1;
  _numQuadPs = 2*_spatialOrder-1;

  _w = alloc_1d<REAL>(_numQuadPs);
  _x = alloc_1d<REAL>(_numQuadPs);

  if (useGaussian)
  { // weights and abscissa for Gaussian quadrature
    gauleg<REAL>(-1.0, 1.0, _x-1, _w-1, _numQuadPs);
  }
  else
  { // weights and abscissa for uniform quadrature: this ensures that
    // quadrature points are uniformly placed
    REAL d = 2.0/_spatialOrder;
    _x[0] = -1.0 + 0.5*d;
    for (unsigned i=1; i<_spatialOrder; ++i)
      _x[i] = _x[i-1] + d;
  }

  // allocate memory for legendre polynomials and their derivatives
  _legpol = alloc_2d_c<REAL>(_numCoeffPerElem, _numQuadPs);
  _dlegpol = alloc_2d_c<REAL>(_numCoeffPerElem, _numQuadPs);
  // compute Legendre polynomials and their derivatives at the
  // abscissa locations
  for(unsigned m=0; m<_numCoeffPerElem; m++)
    for(unsigned cc=0; cc<_numQuadPs; cc++)
    {
      _legpol[m][cc] = legendre_p<REAL>(m,_x[cc]);
      _dlegpol[m][cc] = legendre_p_d<REAL>(m,_x[cc]);
    }

  // allocate memory for normalization constants for basis functions
  _Cconst = alloc_1d<REAL>(_numCoeffPerElem);

  // set normalization constants for basis function
  for(unsigned lm=0; lm<_numCoeffPerElem; ++lm)
    _Cconst[lm] = 1./(2.0*lm+1);

  // allocate memory for normalization constants for basis functions in 2D
  _Cconst2D = alloc_2d_c<REAL>(_numCoeffPerElem, _numCoeffPerElem);

  // set normalization constants for basis function
  for(unsigned lm=0; lm<_numCoeffPerElem; ++lm)
    for(unsigned ln=0; ln<_numCoeffPerElem; ++ln)
      _Cconst2D[lm][ln] = 1./((2.0*lm+1)*(2.0*ln+1));

  // allocate memory for normalization constants for basis functions in 3D
  _Cconst3D = alloc_3d_c<REAL>(_numCoeffPerElem, _numCoeffPerElem, _numCoeffPerElem);

  // set normalization constants for basis function
  for(unsigned lm=0; lm<_numCoeffPerElem; ++lm)
    for(unsigned ln=0; ln<_numCoeffPerElem; ++ln)
      for(unsigned lp=0; lp<_numCoeffPerElem; ++lp)
        _Cconst3D[lm][ln][lp] = 1./((2.0*lm+1)*(2.0*ln+1)*(2.0*lp+1));

}

template <typename REAL>
WxDGQuadrature<REAL>::~WxDGQuadrature()
{
  delete [] _w;
  delete [] _x;
  delete [] _Cconst;
  free_2d_c(_legpol, _numCoeffPerElem, _numQuadPs);
  free_2d_c(_dlegpol, _numCoeffPerElem, _numQuadPs);
  free_2d_c(_Cconst2D, _spatialOrder, _spatialOrder);
  free_3d_c(_Cconst3D, _spatialOrder, _spatialOrder, _spatialOrder);
}

template <typename REAL>
void
WxDGQuadrature<REAL>::evalExpansionLower1D(const REAL coeffs[], REAL val[])
{
  // zap existing values
  for (unsigned me=0; me<_meqn; ++me)
    val[me] = 0.0;

  double pr = 1.0;
  // accumulate
  for (unsigned r=0; r<_numCoeffPerElem; ++r)
  {
    for (unsigned me=0; me<_meqn; ++me)
      val[me] += pr*coeffs[_ind1D.index(me,r)];
    pr = -1*pr;
  }
}

template <typename REAL>
void
WxDGQuadrature<REAL>::evalExpansionUpper1D(const REAL coeffs[], REAL val[])
{
  // zap existing values
  for (unsigned me=0; me<_meqn; ++me)
    val[me] = 0.0;

  // accumulate
  for (unsigned r=0; r<_numCoeffPerElem; ++r)
    for (unsigned me=0; me<_meqn; ++me)
      val[me] += coeffs[_ind1D.index(me,r)];
}

template <typename REAL>
void
WxDGQuadrature<REAL>::evalExpansion1D(int lm,
  const REAL coeffs[], REAL res[])
{
  // zap existing values
  for (unsigned me=0; me<_meqn; ++me)
    res[me] = 0.0;

  // accumulate
  for (unsigned r=0; r<_numCoeffPerElem; ++r)
    for (unsigned me=0; me<_meqn; ++me)
      res[me] += _legpol[r][lm]*coeffs[_ind1D.index(me,r)];
}

template <typename REAL>
void
WxDGQuadrature<REAL>::evalExpansionLower2D(unsigned dir, unsigned l, const REAL coeffs[], REAL val[])
{
  // zap existing values
  for (unsigned me=0; me<_meqn; ++me)
    val[me] = 0.0;

  double pr, pm;

  if (dir == 0)
  {
    // accumulate
    pr = 1.0;
    for (unsigned r=0; r<_spatialOrder; ++r)
    {
      for (unsigned m=0; m<_spatialOrder; ++m)
        for (unsigned me=0; me<_meqn; ++me)
          val[me] += pr*_legpol[m][l]*coeffs[_ind2D.index(me,r,m)];
      pr = -1.0*pr;
    }
  }
  if (dir == 1)
  {
    // accumulate
    pm = 1.0;
    for (unsigned m=0; m<_spatialOrder; ++m)
    {
      for (unsigned r=0; r<_spatialOrder; ++r)
        for (unsigned me=0; me<_meqn; ++me)
          val[me] += _legpol[r][l]*pm*coeffs[_ind2D.index(me,r,m)];
      pm = -1.0*pm;
    }
  }
}

template <typename REAL>
void
WxDGQuadrature<REAL>::evalExpansionUpper2D(unsigned dir, unsigned l, const REAL coeffs[], REAL val[])
{
  // zap existing values
  for (unsigned me=0; me<_meqn; ++me)
    val[me] = 0.0;

  if (dir == 0)
  {
    // accumulate
    for (unsigned r=0; r<_spatialOrder; ++r)
      for (unsigned m=0; m<_spatialOrder; ++m)
        for (unsigned me=0; me<_meqn; ++me)
          val[me] += _legpol[m][l]*coeffs[_ind2D.index(me,r,m)];
  }
  if (dir == 1)
  {
    // accumulate
    for (unsigned r=0; r<_spatialOrder; ++r)
      for (unsigned m=0; m<_spatialOrder; ++m)
        for (unsigned me=0; me<_meqn; ++me)
          val[me] += _legpol[r][l]*coeffs[_ind2D.index(me,r,m)];
  }
}

template <typename REAL>
void
WxDGQuadrature<REAL>::evalExpansion2D(int lm, int ln,
  const REAL coeffs[], REAL res[])
{
  // zap existing values
  for (unsigned me=0; me<_meqn; ++me)
    res[me] = 0.0;

  // accumulate
  for (unsigned r=0; r<_spatialOrder; ++r)
    for (unsigned m=0; m<_spatialOrder; ++m)
      for (unsigned me=0; me<_meqn; ++me)
        res[me] += _legpol[r][lm]*_legpol[m][ln]*coeffs[_ind2D.index(me,r,m)];
}

template <typename REAL>
void
WxDGQuadrature<REAL>::evalExpansionLower3D(unsigned dir, unsigned l, unsigned s,
  const REAL coeffs[], REAL val[])
{
  // zap existing values
  for (unsigned me=0; me<_meqn; ++me)
    val[me] = 0.0;

  double pm, pn, pp;

  if (dir == 0)
  {
    // accumulate
    pm = 1.0;
    for (unsigned m=0; m<_spatialOrder; ++m)
    {
      for (unsigned n=0; n<_spatialOrder; ++n)
      {
        for (unsigned p=0; p<_spatialOrder; ++p)
          for (unsigned me=0; me<_meqn; ++me)
            val[me] += pm*_legpol[n][l]*_legpol[p][s]*coeffs[_ind3D.index(me,m,n,p)];
      }
      pm = -1.0*pm;
    }
  }
  if (dir == 1)
  {
    // accumulate
    pn = 1.0;
    for (unsigned n=0; n<_spatialOrder; ++n)
    {
      for (unsigned m=0; m<_spatialOrder; ++m)
      {
        for (unsigned p=0; p<_spatialOrder; ++p)
          for (unsigned me=0; me<_meqn; ++me)
            val[me] += _legpol[m][l]*pn*_legpol[p][s]*coeffs[_ind3D.index(me,m,n,p)];
      }
      pn = -1.0*pn;
    }
  }
  if (dir == 2)
  {
    // accumulate
    pp = 1.0;
    for (unsigned p=0; p<_spatialOrder; ++p)
    {
      for (unsigned m=0; m<_spatialOrder; ++m)
      {
        for (unsigned n=0; n<_spatialOrder; ++n)
          for (unsigned me=0; me<_meqn; ++me)
            val[me] += _legpol[m][l]*_legpol[n][s]*pp*coeffs[_ind3D.index(me,m,n,p)];
      }
      pp = -1.0*pp;
    }
  }
}

template <typename REAL>
void
WxDGQuadrature<REAL>::evalExpansionUpper3D(unsigned dir, unsigned l, unsigned s,
  const REAL coeffs[], REAL val[])
{
  // zap existing values
  for (unsigned me=0; me<_meqn; ++me)
    val[me] = 0.0;

  if (dir == 0)
  {
    // accumulate
    for (unsigned m=0; m<_spatialOrder; ++m)
      for (unsigned n=0; n<_spatialOrder; ++n)
        for (unsigned p=0; p<_spatialOrder; ++p)
          for (unsigned me=0; me<_meqn; ++me)
            val[me] += _legpol[n][l]*_legpol[p][s]*coeffs[_ind3D.index(me,m,n,p)];
  }
  if (dir == 1)
  {
    // accumulate
    for (unsigned n=0; n<_spatialOrder; ++n)
      for (unsigned m=0; m<_spatialOrder; ++m)
        for (unsigned p=0; p<_spatialOrder; ++p)
          for (unsigned me=0; me<_meqn; ++me)
            val[me] += _legpol[m][l]*_legpol[p][s]*coeffs[_ind3D.index(me,m,n,p)];
  }
  if (dir == 2)
  {
    // accumulate
    for (unsigned p=0; p<_spatialOrder; ++p)
      for (unsigned m=0; m<_spatialOrder; ++m)
        for (unsigned n=0; n<_spatialOrder; ++n)
          for (unsigned me=0; me<_meqn; ++me)
            val[me] += _legpol[m][l]*_legpol[n][s]*coeffs[_ind3D.index(me,m,n,p)];
  }
}

template <typename REAL>
void
WxDGQuadrature<REAL>::evalExpansion3D(int lm, int ln, int lp,
  const REAL coeffs[], REAL res[])
{
  // zap existing values
  for (unsigned me=0; me<_meqn; ++me)
    res[me] = 0.0;

  // accumulate
  for (unsigned m=0; m<_spatialOrder; ++m)
    for (unsigned n=0; n<_spatialOrder; ++n)
      for (unsigned p=0; p<_spatialOrder; ++p)
        for (unsigned me=0; me<_meqn; ++me)
          res[me] +=
            _legpol[m][lm]*_legpol[n][ln]*_legpol[p][lp]*coeffs[_ind3D.index(me,m,n,p)];
}

// instantiations
template class WxDGQuadrature<float>;
template class WxDGQuadrature<double>;
