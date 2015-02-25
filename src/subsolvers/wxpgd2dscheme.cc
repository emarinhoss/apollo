#include "wxpgd2dscheme.h"

// WarpX lib includes
#include <wxcreator.h>
#include <wxlogger.h>
#include <wxlogstream.h>

# include <wxmpimsg.h>

// std includes
#include <vector>
#include <wxmath.h>
#include <limits>
#include <string>
#include <iostream>
#include <cmath>

template <typename REAL>
WxpDG2Dscheme<REAL>::~WxpDG2Dscheme() {
    delete [] _qM;
    delete [] _qP;
    delete [] _fM;
    delete [] _fP;
    delete [] _gM;
    delete [] _gP;
    delete [] _src;
    delete [] _df;
    delete [] _amdq;
    delete [] _apdq;
    delete [] _sx;
    delete [] _sy;
    delete [] _qauxM;
    delete [] _qauxP;
    delete [] _filterdiag;
    delete [] _filterMatrix;
    free_2d_c(_wave,_meqn,_mwave);
}

template <typename REAL>
void
WxpDG2Dscheme<REAL>::setup(const WxCryptSet& wxc, DM dm)
{
  // call base class setup first
  ApSubSolver<REAL>::setup(wxc, dm);

  _polyOrder = wxc.template get<int>("polynomialOrder");
  _cfl = wxc.template get<REAL>("cfl");

  _dirs[0] = 0;
  _dirs[1] = 1;
  _dirs[2] = 2;
  // check if this simulation is a radial one
  if (wxc.has("isRadial"))
  {
    std::string flg = wxc.template get<std::string>("isRadial");
    if (flg == "true") _dirs[1] = 2;
  }

  // set directions for equation system set
  _eqnSet.setDirs(_dirs);
  _eqnSet.setDim(2); // 2D problem
  _eqnSet.setup(wxc);

  // allocate memory for various arrays
  _meqn = _eqnSet.totalEqns();
  _mwave = _eqnSet.totalWaves();

  _srcSet.setNumEqns(_meqn); // set no of equations
  _srcSet.setup(wxc);

  _qM = alloc_1d<REAL>(_meqn);
  _qP = alloc_1d<REAL>(_meqn);
  _qauxM = alloc_1d<REAL>(_meqn);
  _qauxP = alloc_1d<REAL>(_meqn);
  _df = alloc_1d<REAL>(_meqn); // jump
  _src = alloc_1d<REAL>(_meqn);
  _fM = alloc_1d<REAL>(_meqn);
  _fP = alloc_1d<REAL>(_meqn);
  _gM = alloc_1d<REAL>(_meqn);
  _gP = alloc_1d<REAL>(_meqn);
  // allocate memory for waves, speeds and fluctuations
  _apdq = alloc_1d<REAL>(_meqn); // positive fluctuation
  _amdq = alloc_1d<REAL>(_meqn); // negative fluctuation
  _sx = alloc_1d<REAL>(_mwave); // wave speeds
  _sy = alloc_1d<REAL>(_mwave); // wave speeds
  _wave = alloc_2d_c<REAL>(_meqn, _mwave); // waves

  _dataStruct = wxc.template get<std::vector<WxAny> >("DataStructure");
  // add total number of components
  // _dataStruct.push_back((_polyOrder+1)*(_polyOrder+2)/2*_meqn);

  // create function pointer for initial condition
  const WxCryptSet& initCS = wxc.getSet("InitialCondition");
  std::string kind;
  kind = initCS.template get<std::string>("Kind");
  _initFunc = WxCreatorMap<WxFunction<REAL> >::getNew(kind);
  // setup this function
  _initFunc->setup(initCS);

  _quad = new WxpDGGeometry<REAL>(dm, _meqn, _polyOrder);

  // read list of BC subsolvers
//  std::vector<WxAny> bcs;
//  bcs = wxc.template get<std::vector<WxAny> >("boundaryConditions");
//  std::vector<WxAny>::const_iterator i;
//  for (i=bcs.begin(); i!=bcs.end(); ++i)
//      _bcSubSolvers.push_back( wx_any_cast<std::string>(*i) );

  // read filter information
  _orderSP = wxc.template get<int>("filterOrder");
  _cutoff = wxc.template get<REAL>("filterCutoff");
  int Number = _quad->NpElem()*_quad->NpElem();
  _filterdiag = alloc_1d<REAL>(_quad->NpElem());
  _filterMatrix = alloc_1d<REAL>(Number);
  for(unsigned i=0; i<_quad->NpElem(); i++)
      _filterdiag[i] = 1.0;

  int sk = 0;
  for(unsigned i=0; i<_polyOrder+1;i++)
      for(unsigned j=0;j<_polyOrder-i+1;j++)
      {
          if(i+j>=_polyOrder){
              //_filterdiag[sk] = exp(pow(-36.*((i+j-_cutoff)/(_polyOrder-_cutoff)),_orderSP));
              _filterdiag[sk] = _cutoff;
          }
          sk+=1;
      }

}

template <typename REAL>
void
WxpDG2Dscheme<REAL>::init(PetscReal newDt, Vec out)
{
    DM dm;
    VecGetDM(out, &dm);
    PetscSection stateSection;
    DMGetDefaultSection(dm, &stateSection);

    PetscScalar *x;

    PetscInt kStart, kEnd, k, kEndInterior;

    // Independent variables (t,x,y,z)
    REAL txo[5];
    txo[0] = this->getCurrentTime();

    // results returned by the initialization function
    REAL *d = new REAL[_meqn];

    // ****
    // Initialize the solution vector
    // ****
    // Get cells in this processor
    DMPlexGetHeightStratum(dm, 0, &kStart, &kEnd);
    DMPlexGetHybridBounds(dm, &kEndInterior, NULL, NULL, NULL);
    VecGetArray(out, &x);

    for (k = kStart; k < kEndInterior; ++k)
    {
        for(unsigned node=0; node<_quad->NpElem(); node++){
            txo[1] = _quad->Xcoordinate(k,node);
            txo[2] = _quad->Ycoordinate(k,node);

            _initFunc->func(3, txo, d);
            PetscScalar *xc;

            // reference this cell to the proper location on the solution
            // vector
            DMPlexPointLocalRef(dm,k,x,&xc);
            // assign value returned by the initialization function
            // to the solution vector
            if(xc){
                for(unsigned kk=0; kk<_meqn; kk++)
                    xc[node*_meqn+kk] = d[kk];}
        }
    }
    VecRestoreArray(out, &x);

    isInfinityOrNAN(out, "NAN/INF in initialization");

    // Build filtering Matrix
    _quad->calculateFilter(_filterdiag,_filterMatrix);
}

template <typename REAL>
WxStepperStatus<REAL>
WxpDG2Dscheme<REAL>::step(REAL t, REAL dt, Vec in, Vec out)
{
    isInfinityOrNAN(in, "NAN/INF in input Vector to DG step-function");

    WxStepperStatus<REAL> status;
    DM dm;
    VecGetDM(in, &dm);
    PetscScalar *u;
    PetscScalar *ot, *rhs;
    REAL maxSpeed=0.0;

    // local vectors
    Vec locU, locRHS;
    // create local vector
    DMGetLocalVector(dm, &locU);
    DMGetLocalVector(dm, &locRHS);

    // zero entries of the vectors that will be used to store
    // information
    VecZeroEntries(locU);
    VecZeroEntries(locRHS);
    VecZeroEntries(out);

    // get local values of the global vector in into locX
    DMGlobalToLocalBegin(dm, in, INSERT_VALUES, locU);
    DMGlobalToLocalEnd(dm, in, INSERT_VALUES, locU);

    // get start and end of faces
    PetscInt kStart, kEnd, kEndInterior;
    DMPlexGetHeightStratum(dm, 0, &kStart, &kEnd);
    DMPlexGetHybridBounds(dm, &kEndInterior, NULL , NULL, NULL);
    VecGetArray(locU, &u);
    VecGetArray(out, &ot);

    int NpF = _quad->NpFaces(); // Number of nodes per Face
    int NpE = _quad->NpElem(); // Number of nodes per Element
    int NfE = _quad->NfElem(); // Number of faces per Element
    int f_Fmask[NpF*NfE];
    int nM, nP;
    _quad->returnFmask(f_Fmask);

    for(unsigned k=kStart; k<kEndInterior; k++)
    {
        PetscScalar *qVal, *qOut;
        REAL normals[3*NfE], xc[4], geom[5], Fscale[NfE];
        int connect[2*NfE];
        REAL num_flux[NfE*NpF*_meqn];
        REAL fluxRHS[NpE*_meqn], volumeRHS[NpE*_meqn], Gflux[NpE*_meqn], Fflux[NpE*_meqn], SolQ[NpE*_meqn];
        for(unsigned ke=0; ke<NpE*_meqn; ke++)
            SolQ[ke] = 0.0;

        DMPlexPointLocalRef(dm, k, u, &qVal);
        // Element geometric factors
        _quad->GeometricFactors2d(k,geom);
        // Element face normals
        _quad->Normals2d(k,normals);

        // calculate scaling (face length)/(element Jacobian)
        for(unsigned sc=0; sc<NfE; sc++)
            Fscale[sc] = normals[sc*NfE+2]/geom[4];

        // Element to Element to Faces connection
        _quad->ElementTOElementANDFace(k,connect);

//        int k1 = connect[0];
//        int f1 = connect[1];
//        int k2 = connect[2];
//        int f2 = connect[3];
//        int k3 = connect[4];
//        int f3 = connect[5];

        // =========== Compute n*Flux ===========
        for(unsigned F=0; F<NfE; F++)
        {
            for(unsigned nodes=0; nodes<NpF; nodes++)
            {
                DMPlexPointLocalRef(dm, connect[2*F], u, &qOut);
                int F2 = connect[2*F+1];
                for(unsigned comp=0; comp<_meqn; comp++)
                {
                    // ***** Problem with parallel run is happening here ****
                    // Segmentation Violation, probably memory access out of range
                    // ******************************************************
                        nM  = f_Fmask[F*NpF+nodes];
                        _qM[comp] = qVal[nM*_meqn+comp];
                }

                // Apply Boundary conditions
                if(connect[2*F+1]<0)
                {
                    // Advection
                    _qP[0] = _qM[0];
                    // Maxwell
//                    _qP[0] = 0.0;
//                    _qP[1] = 0.0;
//                    _qP[2] = -_qM[2];
//                    _qP[3] = _qM[3];
//                    _qP[4] = _qM[4];
//                    _qP[5] = _qM[5];

                    // Isentropic
//                    REAL xo = 5.0, yo = 0.0, beta = 5.0, gamma = 1.4;
//                    REAL u = 1., v = 0.;

//                    REAL x = _quad->Xcoordinate(k,nodes);
//                    REAL y = _quad->Ycoordinate(k,nodes);
//                    REAL pi = 3.1416;

//                    REAL xmut = x-u*t, ymvt = y-v*t;
//                    REAL r = sqrt(pow((xmut-xo),2) + pow((ymvt-yo),2));

//                    u   = u - beta*exp(1-pow(r,2))*(ymvt-yo)/(2*pi);
//                    v   = v + beta*exp(1-pow(r,2))*(xmut-xo)/(2*pi);
//                    REAL rho1 = pow(1. - ((gamma-1.)*beta*beta*exp(2.*(1.-r*r))/(16.*gamma*pi*pi)),1./(gamma-1.));
//                    REAL p1   = pow(rho1,gamma);

//                    _qP[0] = rho1;
//                    _qP[1] = rho1*u;
//                    _qP[2] = rho1*v;
//                    _qP[3] = 0.0;
//                    _qP[4] = p1/(gamma-1.) + 0.5*rho1*(u*u+v*v);
                }
                else
                {
                    for(unsigned comp=0; comp<_meqn; comp++)
                    {
                        nP = f_Fmask[F2*NpF+NpF-1-nodes];
                        _qP[comp] = qOut[nP*_meqn+comp];
                    }

                }

                // evaluate fluxes
                _eqnSet.flux(0, xc, _qM, _qauxM, _fM);
                _eqnSet.flux(0, xc, _qP, _qauxP, _fP);
                _eqnSet.flux(1, xc, _qM, _qauxM, _gM);
                _eqnSet.flux(1, xc, _qP, _qauxP, _gP);
                // compute jump in Q (q-wave)
                for (unsigned m=0; m<_meqn; ++m)
                    _df[m] = _qP[m] - _qM[m];

                // call Riemann problem solver to get flucuations
                _eqnSet.riemann(0, xc, xc, _qM, _qP, 0, 0, _df, _wave, _sx, _amdq, _apdq);
                _eqnSet.riemann(1, xc, xc, _qM, _qP, 0, 0, _df, _wave, _sy, _amdq, _apdq);

                // compute the fastest propagating wave speed
                REAL lambda = 0.;
                for (unsigned mw=0; mw<_mwave; ++mw)
                    lambda = dmax(lambda, _sx[mw]*_sx[mw], _sy[mw]*_sy[mw]);
                lambda = sqrt(lambda);

                // find maximum propagation speed in the entire domain
                maxSpeed = dmax(maxSpeed,lambda);

                // Lax-Frederick fluxes
                for(unsigned comp=0; comp<_meqn; comp++)
                    num_flux[(F*NpF+nodes)*_meqn+comp] = 0.5*(normals[NfE*F]*(_fM[comp]+_fP[comp])
                            + normals[NfE*F+1]*(_gM[comp]+_gP[comp]) + lambda*(_qM[comp]-_qP[comp]))*Fscale[F];
            }
        }

        // LIFT Fluxes
        _quad->LIFT_flux(k,fluxRHS,num_flux);

        // =========== Compute Volume Integrals ===========
        for(unsigned nodes=0; nodes<NpE; nodes++)
        {
            for(unsigned comp=0; comp<_meqn; comp++)
                _qM[comp] = qVal[nodes*_meqn+comp];

            _eqnSet.flux(0, xc, _qM, _qauxM, _fM);
            _eqnSet.flux(1, xc, _qM, _qauxM, _gM);

            for(unsigned comp=0; comp<_meqn; comp++){
                Gflux[nodes*_meqn+comp] = _gM[comp];
                Fflux[nodes*_meqn+comp] = _fM[comp];}
        }
        // Calculate Weak Derivatives
        _quad->weakDericatives(k,volumeRHS,Fflux,Gflux);

        // filter solution
        for(unsigned i=0; i<NpE;i++)
            for(unsigned j=0; j<NpE; j++)
                for(unsigned me=0; me<_meqn; me++)
                    SolQ[i*NpE+me] += _filterMatrix[i*NpE+j]*(volumeRHS[j*_meqn+me]-fluxRHS[j*_meqn+me]);

        // add all contributions to conserved variable
        DMPlexPointLocalRef(dm,k,ot,&rhs);
        for(unsigned kne=0; kne<NpE*_meqn; kne++){
//            rhs[kne] = SolQ[kne];
            rhs[kne] = volumeRHS[kne]-fluxRHS[kne];
        }
    }

    DMRestoreLocalVector(dm, &locU);
    VecRestoreArray(out, &ot);
    //VecView(out,PETSC_VIEWER_STDOUT_WORLD);
    isInfinityOrNAN(out, "NAN/INF encountered in RHS Vector of DG step-function");

    REAL newDt =  2./3.*_cfl*_quad->dtscale2D()*(_quad->rMin()/maxSpeed);
    status.setStatus(true);
    status.setSuggestedDt(newDt);
    return status;
}

template<typename REAL>
void
WxpDG2Dscheme<REAL>::applyBc(WxpDGGeometry<REAL> quad, REAL dt, Vec inOut)
{
  // apply boundary conditions get to get correct values in ghost
  // cells
  std::vector<std::string>::const_iterator ssi;
  for (ssi=_bcSubSolvers.begin(); ssi!=_bcSubSolvers.end(); ++ssi) {
    // get hold of our sibling BC subsolver
    ApSubSolver<REAL>* ss = this->getParent()->getSubSolver( *ssi );
    // ss->setCurrentTime(Out);
    // cast this to the a grid BC and call step function
    dynamic_cast<WxGridBC<REAL>* >(ss)->applyToArray(quad, dt, inOut);
  }
}

template <typename REAL>
PetscErrorCode
WxpDG2Dscheme<REAL>::isInfinityOrNAN(Vec f, std::string location)
{
    PetscReal fnorm;
    PetscErrorCode ierr;
    ierr = VecNormBegin(f,NORM_2,&fnorm);CHKERRQ(ierr);	/* fnorm <- ||F||  */
    ierr = VecNormEnd(f,NORM_2,&fnorm);CHKERRQ(ierr);
    if (PetscIsInfOrNanReal(fnorm))
    {
        //VecView(f,PETSC_VIEWER_STDOUT_WORLD);
        //REAL test = 0.0;
        WxLogger *l = WxLogger::get("apollo-root.console");
        WxLogStream errStrm = l->getErrorStream();
        errStrm << location ;
        exit(1); // abort execution
    }
    return 0;
}

// instantiations
template class WxpDG2Dscheme<float>;
template class WxpDG2Dscheme<double>;
