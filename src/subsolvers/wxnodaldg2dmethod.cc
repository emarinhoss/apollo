#include "wxnodaldg2dmethod.h"

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
WxNodalDG2dMethod<REAL>::~WxNodalDG2dMethod() {
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
    //delete [] _filterdiag;
    //delete [] _filterMatrix;
    free_2d_c(_wave,_meqn,_mwave);
    delete _geom;
    delete _cub;
    delete _initFunc;
    VecDestroy(&locU);
    VecDestroy(&locRHS);
    DMDestroy(&_dm);
}

template <typename REAL>
void
WxNodalDG2dMethod<REAL>::setup(const WxCryptSet& wxc, DM dm)
{
    _dm = dm;
  // call base class setup first
  ApSubSolver<REAL>::setup(wxc, _dm);

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

  // Calculate connectivity, coordinates and Matrices
  _geom = new wxNodalDGgeometry2D<REAL>(_dm, _meqn, _polyOrder);

  // Evaluate cubatures, surface Gaussian points and derivatives
  _cub = new WxCubature2d<REAL>(_dm,_meqn,_polyOrder,_geom->inverseVandermonde());

  // read list of BC subsolvers
//  std::vector<WxAny> bcs;
//  bcs = wxc.template get<std::vector<WxAny> >("boundaryConditions");
//  std::vector<WxAny>::const_iterator i;
//  for (i=bcs.begin(); i!=bcs.end(); ++i)
//      _bcSubSolvers.push_back( wx_any_cast<std::string>(*i) );

}

template <typename REAL>
void
WxNodalDG2dMethod<REAL>::init(PetscReal newDt, Vec out)
{
//    DM dm;
//    VecGetDM(out, &dm);
    PetscSection stateSection;
    DMGetDefaultSection(_dm, &stateSection);
    REAL maxSpeed = 0., xcc[4];

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
    DMPlexGetHeightStratum(_dm, 0, &kStart, &kEnd);
    DMPlexGetHybridBounds(_dm, &kEndInterior, NULL, NULL, NULL);
    VecGetArray(out, &x);

    for (k = kStart; k < kEnd; ++k)
    {
        for(unsigned node=0; node<_geom->NpElem(); node++){
            txo[1] = _geom->Xcoordinate(k,node);
            txo[2] = _geom->Ycoordinate(k,node);

            _initFunc->func(3, txo, d);
            PetscScalar *xc;

            // reference this cell to the proper location on the solution
            // vector
            DMPlexPointLocalRef(_dm,k,x,&xc);
            // assign value returned by the initialization function
            // to the solution vector
            if(xc){
                for(unsigned kk=0; kk<_meqn; kk++)
                    xc[node*_meqn+kk] = d[kk];}

            // call Riemann problem solver to get flucuations
            _eqnSet.riemann(0, xcc, xcc, d, d, 0, 0, _df, _wave, _sx, _amdq, _apdq);
            _eqnSet.riemann(1, xcc, xcc, d, d, 0, 0, _df, _wave, _sy, _amdq, _apdq);

            // compute the fastest propagating wave speed
            REAL lambda = 0.;
            for (unsigned mw=0; mw<_mwave; ++mw)
                lambda = dmax(lambda, _sx[mw]*_sx[mw], _sy[mw]*_sy[mw]);
            lambda = sqrt(lambda);

            // find maximum propagation speed in the entire domain
            maxSpeed = dmax(maxSpeed,lambda);
        }
    }
    VecRestoreArray(out, &x);

    isInfinityOrNAN(out, "NAN/INF in initialization!");

    // Suggested initial dt
    REAL timeStep = 2./3.*_cfl*_geom->dtscale2D()*(_geom->rMin()/maxSpeed);
    this->setDt(timeStep);

}

template <typename REAL>
WxStepperStatus<REAL>
WxNodalDG2dMethod<REAL>::step(REAL t, REAL dt, Vec in, Vec out)
{
    isInfinityOrNAN(in, "NAN/INF in input Vector to DG step-function");

    WxStepperStatus<REAL> status;
//    DM dm;
//    VecGetDM(in, &dm);
    PetscScalar *u;
    PetscScalar *ot, *rhs;
    REAL maxSpeed=0.0;


    // create local vector
    DMGetLocalVector(_dm, &locU);
    DMGetLocalVector(_dm, &locRHS);

    // zero entries of the vectors that will be used to store
    // information
    VecZeroEntries(locU);
    VecZeroEntries(locRHS);
    VecZeroEntries(out);

    // get local values of the global vector in into locX
    DMGlobalToLocalBegin(_dm, in, INSERT_VALUES, locU);
    DMGlobalToLocalEnd(_dm, in, INSERT_VALUES, locU);

    // get first and last element number
    PetscInt kStart, kEnd, kEndInterior;
    DMPlexGetHeightStratum(_dm, 0, &kStart, &kEnd);
    DMPlexGetHybridBounds(_dm, &kEndInterior, NULL , NULL, NULL);
    VecGetArray(locU, &u);
    VecGetArray(out, &ot);

    int NpE = _geom->NpElem(); // Number of nodes per Element
    int NfE = _geom->NfElem(); // Number of edge/faces per Element
    int Ncubature= _cub->numCubaturePoints(); // Number of cubature points
    int Ngauss = _cub->numGaussianPoints(); // Number of Gaussian points per edge/face

    /**
     * Allocate memory for all vectors
     */
    // Coordinates
    REAL xcoord[NpE], ycoord[NpE];
    REAL xc[4]; // coordinates
    int connect[2*NfE]; // connectivity information element-to-element-to-edge

    // Volume integral
    REAL q_vol[NpE*_meqn], vec_rhs[NpE*_meqn], volInt[NpE*_meqn];
    REAL Iq_vol[Ncubature*_meqn], If_vol[Ncubature*_meqn], Ig_vol[Ncubature*_meqn];

    // Surface integral
    REAL qtemp[NpE*_meqn], q_surf[NpE*_meqn];
    REAL QP[NfE*Ngauss*_meqn], QM[NfE*Ngauss*_meqn], numFlux[NfE*Ngauss*_meqn];
    REAL qgtemp[NfE*Ngauss*_meqn];

    PetscScalar *qVal;
    for(unsigned kelem=kStart; kelem<kEnd; kelem++)
    {
        // get coordinates of all nodes
        for(unsigned nodes=0; nodes<NpE; nodes++)
        {
            xcoord[nodes] = _geom->Xcoordinate(kelem,nodes);
            ycoord[nodes] = _geom->Ycoordinate(kelem,nodes);
        }

        // Get node values for this element
        DMPlexPointLocalRef(_dm, kelem, u, &qVal);
        for(unsigned kk=0; kk<NpE*_meqn; kk++)
            q_vol[kk] = qVal[kk];

        // calculate edge normals and element Jacobian
        REAL geoFacts[5], normals[3*NfE];
        _geom->GeometricFactors2d(kelem,geoFacts); // [drdx, dsdx, drdy, dsdy, J]
        _geom->Normals2d(kelem,normals); // [nx_edge1,ny_edge1,length_edge1, nx_edge2, ny_edge2 ...]

        /** *******************************************************
         *  *******************************************************
         *  Evaluate Volume Integral
         *  *******************************************************
         *  *******************************************************
         */

        // interpolate nodes values into cubature points
        _cub->interpolatedTOCubatures(q_vol,Iq_vol);

        // evaluate the fluxes at each cubature point
        REAL Qvar[_meqn], Qvaraux[_meqn], Fflux[_meqn], Gflux[_meqn];

        for(unsigned point=0; point<Ncubature; point++)
        {

            // Evaluate Flux at each cubature point
            for(unsigned component=0; component<_meqn; component++)
                Qvar[component] = Iq_vol[point*_meqn+component];

            _eqnSet.flux(0, xc, Qvar, Qvaraux, Fflux);
            _eqnSet.flux(1, xc, Qvar, Qvaraux, Gflux);

            for(unsigned component=0; component<_meqn; component++)
            {
                If_vol[point*_meqn+component] = Fflux[component];
                Ig_vol[point*_meqn+component] = Gflux[component];
            }
        }

        // Compute volume terms (dphidx, F) + (dphidy, G)
        _cub->evaluatedVolumeIntegrals(xcoord,ycoord,If_vol,Ig_vol,volInt);

        /** *******************************************************
         *  *******************************************************
         *  Evaluate Surface Integral
         *  *******************************************************
         *  *******************************************************
         */

        // interpolate nodal values to surface Gassian quadrature points
        _cub->nodesTOSurfaceGaussians(q_vol,QM);

        // Flux Gather
        // get the values at the Gaussian points of the adjacent elements
        _geom->ElementTOElementANDFace(kelem,connect);
        for(unsigned edge=0; edge<NfE; edge++)
        {
            // gather only the values for the needed edge
            int edgeNum = connect[2*edge+1];
            // By definition if edgeNum is negative, this is a physical boundary.
            if(edgeNum<0)
            {
                // Apply Boundary conditions
                for(unsigned gpoint=0; gpoint<Ngauss; gpoint++)
                    for(unsigned comp=0; comp<_meqn; comp++)
                        QP[(edge*Ngauss+gpoint)*_meqn+comp] = 0.0;

            }
            else
            {
                // get the values on the element adjacent to this edge
                DMPlexPointLocalRef(_dm, connect[2*edge], u, &qVal);
                for(unsigned kk=0; kk<NpE*_meqn; kk++)
                    qtemp[kk] = qVal[kk];

                // interpolate values at all edges of opposing element
                _cub->nodesTOSurfaceGaussians(qtemp,qgtemp);

                // only use the Gaussian values of the needed edge
                for(unsigned gpoint=0; gpoint<Ngauss; gpoint++)
                    for(unsigned comp=0; comp<_meqn; comp++)
                        QP[(edge*Ngauss+gpoint)*_meqn+comp] = qgtemp[(edgeNum*Ngauss+Ngauss-1-gpoint)*_meqn+comp];
            }
        }

        // Calculate the Numerical Flux
        for(unsigned kk=0; kk<Ngauss*NfE; kk++)
        {
            for(unsigned comp=0; comp<_meqn; comp++)
            {
                _qM[comp] = QM[kk*_meqn+comp];
                _qP[comp] = QP[kk*_meqn+comp];
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
            int curEdge = kk/Ngauss;
            for(unsigned comp=0; comp<_meqn; comp++){
                numFlux[kk*_meqn+comp] = 0.5*(normals[3*curEdge+0]*(_fM[comp]+_fP[comp]) +
                                              normals[3*curEdge+1]*(_gM[comp]+_gP[comp]) +
                                         lambda*(_qM[comp]-_qP[comp]))*normals[3*curEdge+2];
            }
        }

        // Compute surface integral contribution to all nodes
        _cub->calculateSurfaceIntegral(numFlux,q_surf);

        // add surface and volume contributions
        for(unsigned kk=0; kk<NpE*_meqn; kk++)
            volInt[kk] -= q_surf[kk];

        // Multiply by the inverse Mass Matrix
        _geom->multiplyBYinverseMassMatrix(volInt,vec_rhs);

        DMPlexPointLocalRef(_dm,kelem,ot,&rhs);
        for(unsigned kne=0; kne<NpE*_meqn; kne++)
            rhs[kne] = vec_rhs[kne]/geoFacts[4];
    }
    DMRestoreLocalVector(_dm, &locU);
    VecRestoreArray(out, &ot);

    isInfinityOrNAN(out, "NAN/INF encountered in RHS Vector of DG step-function");

    REAL newDt =  2./3.*_cfl*_geom->dtscale2D()*(_geom->rMin()/maxSpeed);
    status.setStatus(true);
    status.setSuggestedDt(newDt);

    return status;
}

template <typename REAL>
PetscErrorCode
WxNodalDG2dMethod<REAL>::isInfinityOrNAN(Vec f, std::string location)
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
template class WxNodalDG2dMethod<float>;
template class WxNodalDG2dMethod<double>;
