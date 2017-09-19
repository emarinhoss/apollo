// WarpX hyper includes
#include "wxhyperboliceqnset.h"

// WarpX lib includes
#include <wxlogger.h>
#include <wxlogstream.h>
#include <wxmath.h>

/**
 * destroy allocated memory
 */
template<typename REAL>
WxHyperbolicEqnSet<REAL>::~WxHyperbolicEqnSet() {
  free_2d_c(_temp_wave,_meqn,_mwave);
  free_2d_c(_temp_rev,_meqn,_meqn);
  free_2d_c(_temp_lev,_meqn,_meqn);
  free_2d_c(_temp_flux,_meqn,_meqn);

  WxHyperbolicEqn<REAL>* tempeqn;
  for (int i=0; i<(int)_eqnSys.size(); ++i)
  {
    tempeqn = _eqnSys[i];
    delete tempeqn;
  }

}

template<typename REAL>
void
WxHyperbolicEqnSet<REAL>::setDirs(unsigned dirs[])
{
  for (unsigned i=0; i<3; ++i)
    _dirs[i] = dirs[i];
}

template<typename REAL>
void
WxHyperbolicEqnSet<REAL>::setDim(unsigned dim)
{
  _ndim = dim;
}

template<typename REAL>
void
WxHyperbolicEqnSet<REAL>::setup(const WxCryptSet& wxc)
{
  // get hold of stream to log debug messages
  WxLogStream infoStrm
    = WxLogger::get("apollo-root.console")->getInfoStream();

  _mwave = _meqn = _mgrads = 0;

  // equation systems to solve
  std::vector<WxAny> eqns = wxc.template get<std::vector<WxAny> >("Equations");
  // loop over each equation adding it to list equations
  std::vector<WxAny>::const_iterator i;
  for (i=eqns.begin(); i!=eqns.end(); ++i)
  {
    // name of cryptset
    std::string eqnName = wx_any_cast<std::string>(*i);
    // find cryptset for equation
    const WxCryptSet& ecs = wxc.getSet(eqnName);
        
    // find its Kind and create equation object
    std::string kind = ecs.template get<std::string>("Kind");
    WxHyperbolicEqn<REAL> *e = WxCreatorMap<WxHyperbolicEqn<REAL> >::getNew(kind);
    // set it up and add it to set of equations
    e->setup(ecs);
    _eqnSys.push_back(e);

    // increment total number of waves and equations
    _mwave += e->mwave();
    _meqn += e->meqn();
    _mgrads += e->mgrads();

    infoStrm << "Equation " <<  eqnName << " is of kind " + kind << std::endl;
  }

  // allocate space for use in Riemann solver
  _temp_wave = alloc_2d_c<REAL>(_meqn, _mwave);

  // allocate space for use in eigenSystem calculation
  _temp_rev = alloc_2d_c<REAL>(_meqn, _meqn);
  _temp_lev = alloc_2d_c<REAL>(_meqn, _meqn);

  // allocate space for flux jabocian
  _temp_flux = alloc_2d_c<REAL>(_meqn, _meqn);

}

template <typename REAL>
void
WxHyperbolicEqnSet<REAL>::rotateToLocalFrame(REAL norm[3], REAL tan1[3], REAL tan2[3], REAL *vin, REAL *vout)
{
  unsigned mloc = 0;
  unsigned meqn;
  typename std::vector<WxHyperbolicEqn<REAL>* >::const_iterator i;
  for (i=_eqnSys.begin(); i!=_eqnSys.end(); ++i)
  {
    meqn = (*i)->meqn();
    // rotate
    (*i)->rotateToLocalFrame(norm, tan1, tan2, vin+mloc, vout+mloc);

    // move location pointer
    mloc += meqn;    
  }  
}

template <typename REAL>
void
WxHyperbolicEqnSet<REAL>::rotateToGlobalFrame(REAL norm[3], REAL tan1[3], REAL tan2[3], REAL *vin, REAL *vout)
{
  unsigned mloc = 0;
  unsigned meqn;
  typename std::vector<WxHyperbolicEqn<REAL>* >::const_iterator i;
  for (i=_eqnSys.begin(); i!=_eqnSys.end(); ++i)
  {
    meqn = (*i)->meqn();
    // rotate
    (*i)->rotateToGlobalFrame(norm, tan1, tan2, vin+mloc, vout+mloc);

    // move location pointer
    mloc += meqn;    
  }  
}

template <typename REAL>
void
WxHyperbolicEqnSet<REAL>::rotateToLocalFrame(REAL norm[3], REAL *vin, REAL *vout)
{
  unsigned mloc = 0;
  unsigned meqn;
  typename std::vector<WxHyperbolicEqn<REAL>* >::const_iterator i;
  for (i=_eqnSys.begin(); i!=_eqnSys.end(); ++i)
  {
    meqn = (*i)->meqn();
    // rotate
    (*i)->rotateToLocalFrame(norm, vin+mloc, vout+mloc);

    // move location pointer
    mloc += meqn;    
  }
}

template <typename REAL>
void
WxHyperbolicEqnSet<REAL>::rotateToGlobalFrame(REAL norm[3], REAL *vin, REAL *vout)
{
  unsigned mloc = 0;
  unsigned meqn;
  typename std::vector<WxHyperbolicEqn<REAL>* >::const_iterator i;
  for (i=_eqnSys.begin(); i!=_eqnSys.end(); ++i)
  {
    meqn = (*i)->meqn();
    // rotate
    (*i)->rotateToGlobalFrame(norm, vin+mloc, vout+mloc);

    // move location pointer
    mloc += meqn;    
  }
}

template<typename REAL>
void
WxHyperbolicEqnSet<REAL>::flux(unsigned d, REAL *x, REAL *q, REAL *qaux, REAL *f)
{
  unsigned mloc = 0;
  unsigned meqn;

  // loop over each equation system, computing fluxes. Fluxes from
  // each equation are accumulated to compute the full flux
  typename std::vector<WxHyperbolicEqn<REAL>* >::const_iterator i;
  for (i=_eqnSys.begin(); i!=_eqnSys.end(); ++i)
  {
    meqn = (*i)->meqn();
    // call flux for the equation
    (*i)->flux(_dirs[d], x, q+mloc, qaux, f+mloc);

    // move location pointer
    mloc += meqn;
  }
}

template<typename REAL>
void
WxHyperbolicEqnSet<REAL>::DGnumericalFlux(REAL *normals, REAL *qM, REAL *qP, REAL *nflux, REAL maxSpeed)
{
  unsigned mloc = 0;
  unsigned meqn;

  // loop over each equation system, computing fluxes. Fluxes from
  // each equation are accumulated to compute the full flux
  typename std::vector<WxHyperbolicEqn<REAL>* >::const_iterator i;
  for (i=_eqnSys.begin(); i!=_eqnSys.end(); ++i)
  {
    meqn = (*i)->meqn();
    // call flux for the equation
    (*i)->DGnumericalFlux(normals, qM+mloc, qP+mloc, nflux+mloc,maxSpeed);

    // move location pointer
    mloc += meqn;
  }
}

template<typename REAL>
void
WxHyperbolicEqnSet<REAL>::DGLimiterTrigger(REAL *qIn, REAL *qOut)
{
  unsigned mloc = 0;
  unsigned meqn;

  // loop over each equation system, computing fluxes. Fluxes from
  // each equation are accumulated to compute the full flux
  typename std::vector<WxHyperbolicEqn<REAL>* >::const_iterator i;
  for (i=_eqnSys.begin(); i!=_eqnSys.end(); ++i)
  {
    meqn = (*i)->meqn();
    // call flux for the equation
    (*i)->DGLimiterTrigger(qIn+mloc, qOut+mloc);

    // move location pointer
    mloc += meqn;
  }
}

template<typename REAL>
void
WxHyperbolicEqnSet<REAL>::fluxJacobian(unsigned d, REAL *x, REAL *q, REAL *qaux, REAL **f)
{
  unsigned mloc = 0;
  unsigned meqn;

  // loop over each equation system, computing fluxes. Fluxes from
  // each equation are accumulated to compute the full flux
  typename std::vector<WxHyperbolicEqn<REAL>* >::const_iterator i;
  for (i=_eqnSys.begin(); i!=_eqnSys.end(); ++i)
  {
    meqn = (*i)->meqn();
    // call flux for the equation
    (*i)->fluxJacobian(_dirs[d], x, q+mloc, qaux, _temp_flux);

    // copy waves into appropriate location
    for (unsigned m=0; m<meqn; ++m)
      for (unsigned mw=0; mw<meqn; ++mw)
        f[m+mloc][mw+mloc] = _temp_flux[m][mw];

    // move location pointer
    mloc += meqn;
  }
}

template<typename REAL>
void
WxHyperbolicEqnSet<REAL>::edgefluxgengeom(unsigned d, REAL *x, REAL *q, REAL *qaux, REAL *f)
{
  unsigned mloc = 0;
  unsigned meqn;

  // loop over each equation system, computing fluxes. Fluxes from
  // each equation are accumulated to compute the full flux
  typename std::vector<WxHyperbolicEqn<REAL>* >::const_iterator i;
  for (i=_eqnSys.begin(); i!=_eqnSys.end(); ++i)
  {
    meqn = (*i)->meqn();
    // call flux for the equation
    (*i)->edgefluxgengeom(_dirs[d], x, q+mloc, qaux, f+mloc);

    // move location pointer
    mloc += meqn;
  }
}

template<typename REAL>
void
WxHyperbolicEqnSet<REAL>::calulateRHS(unsigned N, REAL *geometry, REAL *normals, WxpDGGeometry<REAL> *quad, REAL *q, REAL *dq, REAL *rhs)
{
  unsigned mloc = 0;
  unsigned meqn;

  // loop over each equation system, computing fluxes. Fluxes from
  // each equation are accumulated to compute the full flux
  typename std::vector<WxHyperbolicEqn<REAL>* >::const_iterator i;
  for (i=_eqnSys.begin(); i!=_eqnSys.end(); ++i)
  {
    meqn = (*i)->meqn();
    // call RHS for the equation
    (*i)->RHS(N, geometry, normals, quad, q+mloc, dq+mloc, rhs+mloc);

    // move location pointer
    mloc += meqn*(N+1)*(N+2)/2;
  }
}


template<typename REAL>
void
WxHyperbolicEqnSet<REAL>::riemann(
    unsigned d, REAL *xl, REAL *xr, REAL *ql, REAL *qr, REAL *qauxl, REAL *qauxr,
    REAL *df, REAL **wave, REAL *s, REAL *amdq, REAL *apdq)
{
  unsigned mloc = 0;
  unsigned mwloc = 0;
  unsigned meqn, mwave;

  // loop over each equation system, solving Riemann problem. Waves,
  // speeds and fluctuations from each equation are accumulated to
  // compute the complete Reimann solution
  typename std::vector<WxHyperbolicEqn<REAL>* >::iterator i;
  for (i=_eqnSys.begin(); i!=_eqnSys.end(); ++i)
  {
    meqn = (*i)->meqn(); mwave = (*i)->mwave();
    // set left and right coordinates before calling the Reimann solver
    (*i)->setCoordLeftCell(_ndim, xl);
    (*i)->setCoordRightCell(_ndim, xr);

    // call Riemann solver for the equation
    (*i)->rp(_dirs[d], ql+mloc, qr+mloc, qauxl, qauxr, df+mloc, _temp_wave, s+mwloc, amdq+mloc, apdq+mloc);

    // copy waves into appropriate location
    for (unsigned m=0; m<meqn; ++m)
      for (unsigned mw=0; mw<mwave; ++mw)
        wave[m+mloc][mw+mwloc] = _temp_wave[m][mw];

    // move location pointers
    mloc += meqn;
    mwloc += mwave;
  }
}

template<class REAL>
void
WxHyperbolicEqnSet<REAL>::
riemannt(unsigned td, unsigned d, REAL *xl, REAL *xr, REAL *ql, REAL *qr,
  REAL *amdq, REAL* bmamdq, REAL* bpamdq,
  REAL *apdq, REAL* bmapdq, REAL* bpapdq)
{
  unsigned mloc = 0;
  unsigned meqn;

  // loop over each equation system, solving transverse Riemann
  // problem.
  typename std::vector<WxHyperbolicEqn<REAL>* >::iterator i;
  for (i=_eqnSys.begin(); i!=_eqnSys.end(); ++i)
  {
    meqn = (*i)->meqn();
    // set left and right coordinates before calling the Reimann solver
    (*i)->setCoordLeftCell(_ndim, xl);
    (*i)->setCoordRightCell(_ndim, xr);

    // call transverse Riemann solver for equation
    (*i)->rpt(_dirs[td], _dirs[d], ql+mloc, qr+mloc,
      amdq+mloc, bmamdq+mloc, bpamdq+mloc,
      apdq+mloc, bmapdq+mloc, bpapdq+mloc);

    // move location pointers
    mloc += meqn;
  }  
}


template<class REAL>
void
WxHyperbolicEqnSet<REAL>::
riemanntc(unsigned td, unsigned d, REAL *xl, REAL *xr, REAL *ql, REAL *qr,
  REAL *s, REAL* bms, REAL* bps)
{
  unsigned mloc = 0;
  unsigned meqn;

  // loop over each equation system, solving transverse Riemann
  // problem.
  typename std::vector<WxHyperbolicEqn<REAL>* >::iterator i;
  for (i=_eqnSys.begin(); i!=_eqnSys.end(); ++i)
  {
    meqn = (*i)->meqn();
    // set left and right coordinates before calling the Reimann solver
    (*i)->setCoordLeftCell(_ndim, xl);
    (*i)->setCoordRightCell(_ndim, xr);

    // call transverse Riemann solver for equation
    (*i)->rptc(_dirs[td], _dirs[d], ql+mloc, qr+mloc,
      s+mloc, bms+mloc, bps+mloc);

    // move location pointers
    mloc += meqn;
  }
}


template<typename REAL>
void
WxHyperbolicEqnSet<REAL>::riemannlimit(
    unsigned d, REAL *xl, REAL *xr, REAL *ql, REAL *qr, REAL *ql1, REAL *qr1, REAL *qauxl, 
    REAL *qauxr, REAL *df, REAL **wave, REAL *s, REAL *amdq, REAL *apdq)
{
  unsigned mloc = 0;
  unsigned mwloc = 0;
  unsigned meqn, mwave;

  // loop over each equation system, solving Riemann problem. Waves,
  // speeds and fluctuations from each equation are accumulated to
  // compute the complete Reimann solution
  typename std::vector<WxHyperbolicEqn<REAL>* >::iterator i;
  for (i=_eqnSys.begin(); i!=_eqnSys.end(); ++i)
  {
    meqn = (*i)->meqn(); mwave = (*i)->mwave();
    // set left and right coordinates before calling the Reimann solver
    (*i)->setCoordLeftCell(_ndim, xl);
    (*i)->setCoordRightCell(_ndim, xr);

    // call Riemann solver for the equation
    (*i)->rplimit(_dirs[d], ql+mloc, qr+mloc, ql1+mloc, qr1+mloc, qauxl, qauxr, df+mloc, _temp_wave, s+mwloc, amdq+mloc, apdq+mloc);

    // copy waves into appropriate location
    for (unsigned m=0; m<meqn; ++m)
      for (unsigned mw=0; mw<mwave; ++mw)
        wave[m+mloc][mw+mwloc] = _temp_wave[m][mw];

    // move location pointers
    mloc += meqn;
    mwloc += mwave;
  }
}

template<class REAL>
void
WxHyperbolicEqnSet<REAL>::
eigenSystem(unsigned d, REAL *q, REAL *ev, REAL **lev, REAL **rev)
{
  unsigned mloc = 0;
  unsigned meqn;

  // loop over each equation system, computing eigenvectors and
  // eigenvalues.
  typename std::vector<WxHyperbolicEqn<REAL>* >::iterator i;
  for (i=_eqnSys.begin(); i!=_eqnSys.end(); ++i)
  {
    meqn = (*i)->meqn();

    // call eigenSystem calculator for equation
    (*i)->eigenSystem(_dirs[d], q+mloc, ev+mloc, _temp_lev, _temp_rev);

    // copy left eigenvectors into appropriate location
    for (unsigned m=0; m<meqn; ++m)
      for (unsigned mw=0; mw<meqn; ++mw)
        lev[m+mloc][mw+mloc] = _temp_lev[m][mw];

    // copy right eigenvectors into appropriate location
    for (unsigned m=0; m<meqn; ++m)
      for (unsigned mw=0; mw<meqn; ++mw)
        rev[m+mloc][mw+mloc] = _temp_rev[m][mw];

    // move location pointer
    mloc += meqn;
  }
}

template <typename REAL>
void
WxHyperbolicEqnSet<REAL>::
setFaceVectors(REAL norm[3], REAL tan1[3], REAL tan2[3])
{
  // loop over each equation system, solving Riemann problem. Waves,
  // speeds and fluctuations from each equation are accumulated to
  // compute the complete Reimann solution
  typename std::vector<WxHyperbolicEqn<REAL>* >::iterator i;
  for (i=_eqnSys.begin(); i!=_eqnSys.end(); ++i)
  {
    // call set face vector method for the equation
    (*i)->setFaceVectors(norm, tan1, tan2);
  }
}

template <typename REAL>
void
WxHyperbolicEqnSet<REAL>::
setIndices(int idx[3])
{
  // loop over each equation system, setting the indexing for
  // possible access later in the specific equation system
  typename std::vector<WxHyperbolicEqn<REAL>* >::iterator i;
  for (i=_eqnSys.begin(); i!=_eqnSys.end(); ++i)
  {
    // call set face vector method for the equation
    (*i)->setIndices(idx);
  }
}

template <typename REAL>
void
WxHyperbolicEqnSet<REAL>::
computeConservedAndPrimitiveAVEVariables(int kNodes, PetscScalar *qIn, REAL *AVE, REAL *qCons, REAL *qPrim)
{
    unsigned mloc = 0;
    unsigned meqn;

    // zero entry values for summation
    for(unsigned k=0; k<_meqn; k++)
        qCons[k] = 0.;

    // calculate the average conserved variables
    for(unsigned variables=0; variables<_meqn; variables++)
        for(unsigned nodes=0; nodes<kNodes; nodes++)
            qCons[variables] += AVE[nodes]*qIn[nodes*_meqn+variables];

    // loop over each equation system, computing fluxes. Fluxes from
    // each equation are accumulated to compute the full flux
    typename std::vector<WxHyperbolicEqn<REAL>* >::const_iterator i;
    for (i=_eqnSys.begin(); i!=_eqnSys.end(); ++i)
    {
      meqn = (*i)->meqn();
      // call flux for the equation
      (*i)->primitiveVariables(qCons+mloc,qPrim+mloc);

      // move location pointer
      mloc += meqn;
    }
}

template<typename REAL>
void
WxHyperbolicEqnSet<REAL>::getPrimitiveVariable(REAL *qCons, REAL *qPrim)
{
  unsigned mloc = 0;
  unsigned meqn;

  // loop over each equation system, computing primitives. Primitives from
  // each equation are accumulated to compute the full vector of primitives.
  typename std::vector<WxHyperbolicEqn<REAL>* >::const_iterator i;
  for (i=_eqnSys.begin(); i!=_eqnSys.end(); ++i)
  {
    meqn = (*i)->meqn();
    // call flux for the equation
    (*i)->primitiveVariables(qCons+mloc, qPrim+mloc);

    // move location pointer
    mloc += meqn;
  }
}

template<typename REAL>
void
WxHyperbolicEqnSet<REAL>::tuAndAliabadiLimiter(REAL *avgCons, REAL *avgPrim, REAL *dGrads, REAL *limitedValues)
{
  unsigned mloc = 0;
  unsigned meqn;

  // loop over each equation system, computing fluxes. Fluxes from
  // each equation are accumulated to compute the full flux
  typename std::vector<WxHyperbolicEqn<REAL>* >::const_iterator i;
  for (i=_eqnSys.begin(); i!=_eqnSys.end(); ++i)
  {
    meqn = (*i)->meqn();
    // call flux for the equation
    (*i)->limiterTuAndAliabadi(avgCons+mloc, avgPrim+mloc, dGrads+mloc, limitedValues+mloc);

    // move location pointer
    mloc += meqn;
  }
}

// instantiations
template class WxHyperbolicEqnSet<float>;
template class WxHyperbolicEqnSet<double>;
