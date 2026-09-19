// WarpX hyper include
#include "wx3dfemgeom.h"

// WarpX lib includes
#include <wxmath.h>

template <typename REAL>
Wx3dFEMGeom<REAL>::Wx3dFEMGeom(unsigned meqn, unsigned femSpor, unsigned dgSpor, WxGridBox<REAL> gb)
    : _meqn(meqn), _femSpOr(femSpor), _dgSpOr(dgSpor),
      _ind3D(WxRange(0,_meqn, 0,_dgSpOr, 0,_dgSpOr, 0,_dgSpOr))
{

  _ndims = gb.ndims();
  _nodesPerElem = pow(_femSpOr+1,_ndims);
  _modesPerElem = pow(_dgSpOr,_ndims);
  _quadPoints = ceil((_femSpOr+1)/2+1);

  _nxNodesElem = gb.ncells(0);
  _nyNodesElem = gb.ncells(1);
  _nzNodesElem = gb.ncells(2);
  _TotNumNodes = _nxNodesElem*_nyNodesElem*_nzNodesElem;

  _TotNumElemNX = gb.ncells(0)-1/_femSpOr;
  _TotNumElemNY = gb.ncells(1)-1/_femSpOr;
  _TotNumElemNZ = gb.ncells(2)-1/_femSpOr;
  _TotNumElem  = _TotNumElemNX*_TotNumElemNY*_TotNumElemNX;

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

  getNodeCoordinates(gb, _coords);
  buildConnectivity( _connect);

  // Mapping between ijk coordinates and element numbering
  _ijktoen = alloc_3d_c<int>(_TotNumElemNX,_TotNumElemNY,_TotNumElemNZ);
  _entoijk = alloc_2d_c<int>(_TotNumElem,3);
  ijktoelementnumber(_ijktoen,_entoijk);

  // allocate memory for normalization constants for basis functions in 3D
  _Cconst3D = alloc_3d_c<REAL>(_modesPerElem, _modesPerElem, _modesPerElem);
  // set normalization constants for basis function
  for(unsigned lm=0; lm<_modesPerElem; ++lm)
    for(unsigned ln=0; ln<_modesPerElem; ++ln)
      for(unsigned lp=0; lp<_modesPerElem; ++lp)
        _Cconst3D[lm][ln][lp] = 1./((2.0*lm+1)*(2.0*ln+1)*(2.0*lp+1));
}

template <typename REAL>
Wx3dFEMGeom<REAL>::~Wx3dFEMGeom()
{
  delete [] _w;
  delete [] _x;
  free_2d_c(_coords, _TotNumNodes, _ndims);
  free_2d_c(_connect, _TotNumElem, _nodesPerElem);
  free_2d_c(_legpol, _modesPerElem, _quadPoints);
  free_2d_c(_dlegpol, _modesPerElem, _quadPoints);
  free_2d_c(_Phi, _nodesPerElem, _quadPoints);
  free_2d_c(_dPhi, _nodesPerElem, _quadPoints);
  free_2d_c(_entoijk, _TotNumElem, 3);
  free_3d_c(_Cconst3D, _modesPerElem, _modesPerElem, _modesPerElem);
  free_3d_c(_ijktoen, _TotNumElemNX, _TotNumElemNY, _TotNumElemNZ);

}

template <typename REAL>
void
Wx3dFEMGeom<REAL>::evalQvalue3D(int li, int lj, int lk, REAL **q1, REAL *res)
{
  // zap existing values
  for (unsigned me=0; me<_meqn; ++me)
    res[me] = 0.0;

  // accumulate
  for (unsigned me=0; me<_meqn; ++me)
      for (unsigned nr = 0; nr<_nodesPerElem; ++nr)
          res[me] += _Phi[nr][li]*_Phi[nr][lj]*_Phi[nr][lk]*q1[me][nr];
}


template <typename REAL>
void
Wx3dFEMGeom<REAL>::buildConnectivity(int **rhs)
{
    int NPD = _femSpOr+1; // number of nodes per direction

    // loop over elements
    for(unsigned nz=0; nz<_TotNumElemNZ; nz++)
        for(unsigned ny=0; ny<_TotNumElemNY; ny++)
            for(unsigned nx=0; nx<_TotNumElemNX; nx++)
            {
                // loop over nodes on each element
                for(unsigned k=0; k<NPD; k++)
                    for(unsigned j=0; j<NPD; j++)
                        for(unsigned i=0; i<NPD; i++)
                            rhs[nz*_TotNumElemNZ+ny*_TotNumElemNY+nx][k*NPD*NPD+j*NPD+i] =
                                    (nx*_femSpOr)+(ny*_nxNodesElem*_femSpOr)+(nz*_nxNodesElem*_nyNodesElem*_femSpOr)+
                                    k*_nxNodesElem*_nyNodesElem+j*_nxNodesElem+i;
            }
}

template <typename REAL>
void
Wx3dFEMGeom<REAL>::getNodeCoordinates(WxGridBox<REAL> gb, REAL **rhs)
{
    REAL xc[4];
    unsigned n = 0;

    for(unsigned k=0;k<gb.upperIdx(2);k++)
        for(unsigned j=0; j<gb.upperIdx(1); j++)
            for(unsigned i=0; i<gb.upperIdx(0); i++)
            {
                gb.nodeCoord(i,j,k,xc);
                for(int j=0;j<_ndims;j++)
                {
                    rhs[n][j] = xc[j];
                    n += 1;
                }
            }
}

template <typename REAL>
void
Wx3dFEMGeom<REAL>::ijktoelementnumber(int ***rhs, int **rhs2)
{
    int num = 0;

    for(unsigned k=0; k<_TotNumElemNZ; k++)
        for(unsigned j=0; j<_TotNumElemNY; j++)
            for(unsigned i=0; i<_TotNumElemNX; i++)
            {
                rhs[i][j][k] = num;
                rhs2[num][0] = i;
                rhs2[num][1] = j;
                rhs2[num][2] = k;
                num += 1;
            }
}

// instantiations
//template class Wx3dFEMGeom<float>;
template class Wx3dFEMGeom<double>;
