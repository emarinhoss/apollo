// WarpX hyper include
#include "wxfemgeometry.h"

// WarpX lib includes
#include <wxmath.h>

template <typename REAL>
WxFEMGeometry<REAL>::WxFEMGeometry(unsigned ieqn, unsigned meqn, unsigned femSpor, unsigned dgSpor, WxGridBox<REAL> gb, unsigned nodes)
    : _meqn(meqn), _femSpOr(femSpor), _dgSpOr(dgSpor), _TotNumNodes(nodes), _mauxeqn(ieqn),
      _ind1D(WxRange(0,_mauxeqn, 0,(_dgSpOr))),
      _ind2D(WxRange(0,_mauxeqn, 0,(_dgSpOr), 0,(_dgSpOr))),
      _ind3D(WxRange(0,_mauxeqn, 0,(_dgSpOr), 0,(_dgSpOr), 0,(_dgSpOr)))
{

  _ndims = gb.ndims();
  _nodesPerElem = pow(_femSpOr+1,_ndims);
  _modesPerElem = _dgSpOr*_ndims;
  unsigned spatialO = ( _femSpOr>_dgSpOr ) ? _femSpOr : _dgSpOr;
  _quadPoints = ceil((3*spatialO+1)/2+1);
  _TotNumElem  = (_TotNumNodes-1)/_femSpOr;

  // allocate memory for the Gauss-Legendre quadrature
  _w = alloc_1d<REAL>(_quadPoints); // weights
  _x = alloc_1d<REAL>(_quadPoints); // coordinates

  // calculate weights and abscissa for Gauss-Legendre quadrature
  gauleg<REAL>(-1.0, 1.0, _x-1, _w-1, _quadPoints);

    // allocate memory for the Nodes Global Coordinate and Connectivity matrix
  _coords  = alloc_2d_c<REAL>(_TotNumNodes,_ndims);
  _connect = alloc_2d_c<int>(_TotNumElem, _nodesPerElem);

  // allocate memory for legendre polynomials and their derivatives
  _legpol = alloc_2d_c<REAL>(_modesPerElem, _quadPoints);
  _dlegpol = alloc_2d_c<REAL>(_modesPerElem, _quadPoints);

  // allocate memory for shape functions and their derivatives
  _Phi = alloc_2d_c<REAL>(_nodesPerElem, _quadPoints);
  _dPhi= alloc_2d_c<REAL>(_nodesPerElem, _quadPoints);

  // compute Legendre/Lagrange polynomials and their derivatives at the
  // quadrature locations
  for(unsigned m=0; m<_modesPerElem; m++)
    for(unsigned cc=0; cc<_quadPoints; cc++)
    {
        _legpol[m][cc] = legendre_p<REAL>(m,_x[cc]);
        _dlegpol[m][cc] = legendre_p_d<REAL>(m,_x[cc]);
    }

  for(unsigned m=0; m<_nodesPerElem; m++)
    for(unsigned cc=0; cc<_quadPoints; cc++)
    {
        _Phi[m][cc] = lagrange_p(_femSpOr,_x[cc],m);
        _dPhi[m][cc] = lagrange_p_d(_femSpOr,_x[cc],m);
    }

  //if(_ndims>1){_TotNumElem = _TotNumElem*((gb.ncells(1)-1)/_femSpOr);}
  //if(_ndims>2){_TotNumElem = _TotNumElem*((gb.ncells(2)-1)/_femSpOr);}


  switch (_ndims)
  {
  case 1:
      getCoordinate1d(gb, _coords);
      buildConnectivity1d(_nodesPerElem, _connect);
      break;

//  case 2:
//      getCoordinate2d(gb, _coords);
//      buildConnectivity2d(_ndims, _TotNumNodes, _connect);
//      break;

//  case 3:
//      getCoordinate3d(gb, _coords);
//      buildConnectivity3d(_ndims, _TotNumNodes, _connect);
//      break;

  default:
      REAL test = 0.0;

  }

  // allocate memory for normalization constants for basis functions
  _Cconst = alloc_1d<REAL>(_modesPerElem);
  // set normalization constants for basis function
  for(unsigned lm=0; lm<_modesPerElem; ++lm)
    _Cconst[lm] = 1./(2.0*lm+1);

  // allocate memory for normalization constants for basis functions in 2D
  _Cconst2D = alloc_2d_c<REAL>(_modesPerElem, _modesPerElem);
  // set normalization constants for basis function
  for(unsigned lm=0; lm<_modesPerElem; ++lm)
    for(unsigned ln=0; ln<_modesPerElem; ++ln)
      _Cconst2D[lm][ln] = 1./((2.0*lm+1)*(2.0*ln+1));

  // allocate memory for normalization constants for basis functions in 3D
  _Cconst3D = alloc_3d_c<REAL>(_modesPerElem, _modesPerElem, _modesPerElem);
  // set normalization constants for basis function
  for(unsigned lm=0; lm<_modesPerElem; ++lm)
    for(unsigned ln=0; ln<_modesPerElem; ++ln)
      for(unsigned lp=0; lp<_modesPerElem; ++lp)
        _Cconst3D[lm][ln][lp] = 1./((2.0*lm+1)*(2.0*ln+1)*(2.0*lp+1));
}

template <typename REAL>
WxFEMGeometry<REAL>::~WxFEMGeometry()
{
  delete [] _w;
  delete [] _x;
  delete [] _Cconst;
  free_2d_c(_coords, _TotNumNodes, _ndims);
  free_2d_c(_connect, _TotNumElem, _nodesPerElem);
  free_2d_c(_legpol, _modesPerElem, _quadPoints);
  free_2d_c(_dlegpol, _modesPerElem, _quadPoints);
  free_2d_c(_Phi, _nodesPerElem, _quadPoints);
  free_2d_c(_dPhi, _nodesPerElem, _quadPoints);
  free_2d_c(_Cconst2D, _modesPerElem, _modesPerElem);
  free_3d_c(_Cconst3D, _modesPerElem, _modesPerElem, _modesPerElem);

}

template <typename REAL>
void
WxFEMGeometry<REAL>::evalExpansion1Dcg(int lm, REAL **q1, REAL *res)
{
    // zap existing values
    for (unsigned me=0; me<_meqn; ++me)
      res[me] = 0.0;

    // accumulate
    for (unsigned me=0; me<_meqn; ++me)
        for (unsigned nr = 0; nr<_nodesPerElem; ++nr)
            res[me] += _Phi[nr][lm]*q1[me][nr];
}

template <typename REAL>
void
WxFEMGeometry<REAL>::evalExpansionLower1Ddg(const REAL coeffs[], REAL val[])
{
  // zap existing values
  for (unsigned me=0; me<_mauxeqn; ++me)
    val[me] = 0.0;

  double pr = 1.0;
  // accumulate
  for (unsigned r=0; r<_modesPerElem; ++r)
  {
    for (unsigned me=0; me<_mauxeqn; ++me)
      val[me] += pr*coeffs[_ind1D.index(me,r)];
    pr = -1*pr;
  }
}

template <typename REAL>
void
WxFEMGeometry<REAL>::evalExpansionUpper1Ddg(const REAL coeffs[], REAL val[])
{
  // zap existing values
  for (unsigned me=0; me<_mauxeqn; ++me)
    val[me] = 0.0;

  // accumulate
  for (unsigned r=0; r<_modesPerElem; ++r)
    for (unsigned me=0; me<_mauxeqn; ++me)
      val[me] += coeffs[_ind1D.index(me,r)];
}

template <typename REAL>
void
WxFEMGeometry<REAL>::evalExpansion1Ddg(int lm, const REAL coeffs[], REAL res[])
{
    // zap existing values
    for (unsigned me=0; me<_mauxeqn; ++me)
      res[me] = 0.0;

    // accumulate
    for (unsigned r=0; r<_modesPerElem; ++r)
      for (unsigned me=0; me<_mauxeqn; ++me)
        res[me] += _legpol[r][lm]*coeffs[_ind1D.index(me,r)];
}

template <typename REAL>
void
WxFEMGeometry<REAL>::buildConnectivity1d(int spOrd, int **rhs)
{
    for(int k=0; k<_TotNumElem; k++)
        for(int j=0; j<_nodesPerElem; j++)
            rhs[k][j] = (spOrd-1)*k+j-_femSpOr;
}

template <typename REAL>
void
WxFEMGeometry<REAL>::getCoordinate1d(WxGridBox<REAL> gb, REAL **rhs)
{
    REAL xc[5];
    _dx = (gb.upper(0)-gb.lower(0))/(_TotNumElem-2);

    for(int i=0; i<_TotNumNodes; i++)
    {
        gb.femNodeCoord(i-_femSpOr,xc+1);
        for(int j=0;j<_ndims;j++)
            rhs[i][j] = xc[j+1];
    }
}

// instantiations
template class WxFEMGeometry<float>;
template class WxFEMGeometry<double>;
