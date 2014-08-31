// WarpX hyper include
#include "wxfemquadrature.h"

// WarpX lib includes
#include <wxmath.h>
#include <wxfemshapefuncs.h>

template <typename REAL>
WxFEMQuadrature<REAL>::WxFEMQuadrature(unsigned meqn, unsigned nodesPerElement, WxGridBox<REAL> gb, unsigned qvals)
    : _meqn(meqn), _nodesPerElem(nodesPerElement), _numNodes(qvals),
    _ind1D(WxRange(0,meqn))
{

  _numQuadPs = ceil(_nodesPerElem/2+1);
  _w = alloc_1d<REAL>(_nodesPerElem);
  _dgw = alloc_1d<REAL>(_nodesPerElem-1);

  _x = alloc_1d<REAL>(_nodesPerElem);
  _dgx = alloc_1d<REAL>(_nodesPerElem-1);

  _npad = _nodesPerElem-1;
  if(_nodesPerElem==2)
      _npad = 2;

  //{ // weights and abscissa for Gauss-Legendre quadrature
  gauleg<REAL>(-1.0, 1.0, _x-1, _w-1, _numQuadPs);
  gauleg<REAL>(-1.0, 1.0, _dgx-1, _dgw-1, _numQuadPs-1);

  // allocate memory for shape functions
  _legpol = alloc_2d_c<REAL>(_nodesPerElem, _numQuadPs);
  _dlegpol = alloc_2d_c<REAL>(_nodesPerElem, _numQuadPs);

  _dglegpol = alloc_2d_c<REAL>(_nodesPerElem, _numQuadPs-1);
  _dgdlegpol = alloc_2d_c<REAL>(_nodesPerElem, _numQuadPs-1);

  // compute shape functions and their derivatives at the
  // abscissa locations
  //shapes1D(_nodesPerElem, _numQuadPs, _x, _legpol, _dlegpol);
  for(unsigned m=0; m<_nodesPerElem; m++)
    for(unsigned cc=0; cc<_numQuadPs; cc++)
    {
        _legpol[m][cc] = lagrange_p(_nodesPerElem-1,_x[cc],m);
        _dlegpol[m][cc] = lagrange_p_d(_nodesPerElem-1,_x[cc],m);
    }

  shapes1D(_nodesPerElem, _numQuadPs-1, _dgx, _dglegpol, _dgdlegpol);

  _numElem  = (_numNodes-1)/(_nodesPerElem-1);
  _coords = alloc_1d<REAL>(_numNodes);
  _connect = alloc_2d_c<int>(_numElem, _nodesPerElem);
  getCoordinate(gb, _coords);
  buildConnectivity(_nodesPerElem, _numNodes, _connect);

}

template <typename REAL>
WxFEMQuadrature<REAL>::~WxFEMQuadrature()
{
  delete [] _w;
  delete [] _x;
  delete [] _dgw;
  delete [] _dgx;
  delete [] _coords;
  free_2d_c(_connect, _numElem, _nodesPerElem);
  free_2d_c(_legpol, _nodesPerElem, _numQuadPs);
  free_2d_c(_dlegpol, _nodesPerElem, _numQuadPs);
  free_2d_c(_dglegpol, _nodesPerElem, _numQuadPs-1);
  free_2d_c(_dgdlegpol, _nodesPerElem, _numQuadPs-1);
}

template <typename REAL>
void
WxFEMQuadrature<REAL>::evalExpansion1D(int lm, REAL **q1, REAL *res)
{
  // zap existing values
  for (unsigned me=0; me<_meqn; ++me)
    res[me] = 0.0;

  // accumulate
  for (unsigned me=0; me<_meqn; ++me)
      for (unsigned nr = 0; nr<_nodesPerElem; ++nr)
          res[me] += _legpol[nr][lm]*q1[me][nr];
}

template <typename REAL>
void
WxFEMQuadrature<REAL>::evalExpansion1Ddg(int auxMeqn, int lm, REAL **q1, REAL *res)
{
  // zap existing values
  for (unsigned me=0; me<auxMeqn; ++me)
    res[me] = 0.0;

  // accumulate
  for (unsigned me=0; me<auxMeqn; ++me)
      for (unsigned nr = 0; nr<_nodesPerElem; ++nr)
          res[me] += _dglegpol[nr][lm]*q1[me][nr];
}

template <typename REAL>
void
WxFEMQuadrature<REAL>::buildConnectivity(unsigned spOrd, unsigned n, int **rhs)
{
    unsigned elems = (n-1)/(spOrd-1); //number of FEM elements

    for(unsigned k=0; k<elems; ++k)
        for(unsigned j=0; j<spOrd; ++j){
            rhs[k][j] = (spOrd-1)*k+j-_npad;}
}

template <typename REAL>
void
WxFEMQuadrature<REAL>::getCoordinate(WxGridBox<REAL> gb, REAL *rhs)
{
    REAL xd[5];
    REAL lb = gb.lower(0);
    REAL lu = gb.upper(0);
    //REAL dx_node = (lu-lb)/(_numNodes-1);
    _dx = (lu-lb)/(_numElem-2);
    //_dx = gb.dx(0)*(_nodesPerElem-1);

    for(int i=0; i<_numNodes; ++i){
        //gb.coord(i-_npad, xd);
        gb.femNodeCoord(i-_npad,xd+1);
        //rhs[i] = lb+(i-_npad)*dx_node;
        rhs[i] = xd[1];
    }
}

// instantiations
template class WxFEMQuadrature<float>;
template class WxFEMQuadrature<double>;
