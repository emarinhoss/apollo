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
    delete [] _filterdiag;
    delete [] _filterMatrix;
    free_2d_c(_wave,_meqn,_mwave);
}

template <typename REAL>
void
WxNodalDG2dMethod<REAL>::setup(const WxCryptSet& wxc, DM dm)
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

  _geom = new wxNodalDGgeometry2D<REAL>(dm, _meqn, _polyOrder);

  _cub = new WxCubature2d<REAL>(dm,_meqn,_polyOrder,_geom->inverseVandermonde());

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
    DM dm;
    VecGetDM(out, &dm);
    PetscSection stateSection;
    DMGetDefaultSection(dm, &stateSection);
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
    DMPlexGetHeightStratum(dm, 0, &kStart, &kEnd);
    DMPlexGetHybridBounds(dm, &kEndInterior, NULL, NULL, NULL);
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
            DMPlexPointLocalRef(dm,k,x,&xc);
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

    isInfinityOrNAN(out, "NAN/INF in initialization");

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

    // get first and last element number
    PetscInt kStart, kEnd, kEndInterior;
    DMPlexGetHeightStratum(dm, 0, &kStart, &kEnd);
    DMPlexGetHybridBounds(dm, &kEndInterior, NULL , NULL, NULL);
    VecGetArray(locU, &u);
    VecGetArray(out, &ot);

    int NpF = _geom->NpFaces(); // Number of nodes per Face
    int NpE = _geom->NpElem(); // Number of nodes per Element
    int NfE = _geom->NfElem(); // Number of faces per Element
    int Ncubature= _cub->numCubaturePoints(); // Number of cubature points
    int Ngauss = _cub->numGaussianPoints(); // Number of Gaussian points per edge/face
    int f_Fmask[NpF*NfE];
    _geom->returnFmask(f_Fmask);

    /**
     * Allocate memory for all vectors
     */
    // Coordinates
    Vec xcoord, ycoord;
    VecCreateSeq(PETSC_COMM_SELF,NpE,&xcoord);
    VecDuplicate(xcoord,&ycoord);

    // Volume integral
    Vec q_vol, Iq_vol, If_vol, Ig_vol, volInt, vec_rhs;
    VecCreateSeq(PETSC_COMM_SELF,NpE*_meqn,&q_vol);
    VecDuplicate(q_vol,&vec_rhs);
    VecCreateSeq(PETSC_COMM_SELF,Ncubature*_meqn,&Iq_vol);
    VecDuplicate(Iq_vol,&If_vol);
    VecDuplicate(Iq_vol,&Ig_vol);
    VecDuplicate(Iq_vol,&volInt);

    // Surface integral
    Vec QP, QM, numFlux, qtemp, qgtemp, Flux, q_surf;
    VecCreateSeq(PETSC_COMM_SELF,NpE*_meqn,&qtemp);
    VecCreateSeq(PETSC_COMM_SELF,NfE*Ngauss*_meqn,&qgtemp);
    VecDuplicate(qgtemp,&Flux);
    VecDuplicate(qgtemp,&QP);
    VecDuplicate(qgtemp,&QM);
    VecDuplicate(qgtemp,&numFlux);
    VecDuplicate(qtemp,&q_surf);

    for(unsigned kelem=kStart; kelem<kEnd; kelem++)
    {
        PetscScalar *qVal, *xx, *yy, *zz;
        REAL xc[4];
        int connect[2*NfE];

        // get coordinates of all nodes
        VecGetArray(xcoord,&xx);
        VecGetArray(ycoord,&yy);
        for(unsigned nodes=0; nodes<NpE; nodes++)
        {
            xx[nodes] = _geom->Xcoordinate(kelem,nodes);
            yy[nodes] = _geom->Ycoordinate(kelem,nodes);
        }
        VecRestoreArray(xcoord,&xx);
        VecRestoreArray(ycoord,&yy);

        // Get node values for this element
        DMPlexPointLocalRef(dm, kelem, u, &qVal);
        VecGetArray(q_vol,&zz);
        for(unsigned kk=0; kk<NpE*_meqn; kk++)
            zz[kk] = qVal[kk];
        VecRestoreArray(q_vol,&zz);

        // calculate face normals and element Jacobian
        REAL geoFacts[5], normals[3*NfE], Fscale[NfE];
        _geom->GeometricFactors2d(kelem,geoFacts);
        _geom->Normals2d(kelem,normals);
        for(unsigned kk=0; kk<NfE; kk++)
            Fscale[kk] = normals[kk*NfE+2]/geoFacts[4];


        /**
         * Evaluate Volume Integral
         */

        // interpolate nodes values into cubature points
        _cub->interpolatedTOCubatures(q_vol,Iq_vol);

        // evaluate the fluxes at each cubature point
        VecDuplicate(Iq_vol,&If_vol);
        VecDuplicate(Iq_vol,&Ig_vol);
        REAL Qvar[_meqn], Qvaraux[_meqn], Fflux[_meqn], Gflux[_meqn];

        VecGetArray(Iq_vol,&xx);
        VecGetArray(If_vol,&yy);
        VecGetArray(Ig_vol,&zz);
        for(unsigned kk=0; kk<Ncubature; kk++)
        {

            for(unsigned mm=0; mm<_meqn; mm++)
                Qvar[kk] = xx[kk*_meqn+mm];

            _eqnSet.flux(0, xc, Qvar, Qvaraux, Fflux);
            _eqnSet.flux(1, xc, Qvar, Qvaraux, Gflux);

            for(unsigned mm=0; mm<_meqn; mm++)
            {
                yy[kk*_meqn+mm] = Fflux[mm];
                zz[kk*_meqn+mm] = Gflux[mm];
            }
        }
        VecRestoreArray(Iq_vol,&xx);
        VecRestoreArray(If_vol,&yy);
        VecRestoreArray(Ig_vol,&zz);

        // Compute volume terms (dphidx, F) + (dphidy, G)
        _cub->evaluatedVolumeIntegrals(xcoord,ycoord,If_vol,Ig_vol,&volInt);

        /**
         * Evaluate Surface Integral
         */

        // interpolate nodal values into surface Gassian points
        _cub->nodesTOSurfaceGaussians(q_vol,&QM);

        // Flux Gather
        // get the values at the Gaussian points of the adjacent elements
        _geom->ElementTOElementANDFace(kelem,connect);
        for(unsigned edge=0; edge<NfE; edge++)
        {
            // get the values on the element adjacent to this face
            DMPlexPointLocalRef(dm, connect[2*edge], u, &qVal);
            VecGetArray(qtemp,&xx);
            for(unsigned kk=0; kk<NpE*_meqn; kk++)
                xx[kk] = qVal[kk];
            VecRestoreArray(qtemp,&xx);

            // interpolate values at all faces of opposing element
            _cub->nodesTOSurfaceGaussians(qtemp,&qgtemp);

            // gather only the values for the needed edge
            int edgeNum = connect[2*edge+1];
            VecGetArray(qgtemp,&xx);
            VecGetArray(QP,&yy);
            if(edgeNum<0)
            {
                // Apply Boundary conditions

            }
            else
            {
                for(unsigned gpoint=0; gpoint<Ngauss; gpoint++)
                {
                    for(unsigned comp=0; comp<_meqn; comp++){

                        yy[(edge*Ngauss+gpoint)*_meqn+comp] = xx[(edgeNum*Ngauss+Ngauss-1)*_meqn+comp];
                    }
                }
            }
            VecRestoreArray(qgtemp,&xx);
            VecRestoreArray(QP,&xx);
        }

        // Evaluate the Numerical Flux
        for(unsigned kk=0; kk<Ngauss*NfE; kk++)
        {
            VecGetArray(QM,&xx);
            VecGetArray(QP,&yy);
            for(unsigned comp=0; comp<_meqn; comp++)
            {
                _qM[comp] = xx[kk*_meqn+comp];
                _qP[comp] = yy[kk*_meqn+comp];
            }
            VecRestoreArray(QM,&xx);
            VecRestoreArray(QP,&xx);

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
            VecGetArray(numFlux,&xx);
            for(unsigned comp=0; comp<_meqn; comp++){
                xx[kk*_meqn+comp] = 0.5*(normals[curEdge*NfE]  *(_fM[comp]+_fP[comp]) +
                                         normals[curEdge*NfE+1]*(_gM[comp]+_gP[comp]) +
                                         lambda*(_qM[comp]-_qP[comp]))*Fscale[curEdge];
            }
            VecRestoreArray(numFlux,&xx);
        }

        // Compute surface integral terms
        _cub->calculateSurfaceIntegral(numFlux,q_surf);

        // add surface and volume contributions
        VecAXPY(q_vol,-1.0,q_surf);

        // Multiply by the inverse Mass Matrix
        _geom->multiplyBYinverseMassMatrix(q_vol,vec_rhs);

        DMPlexPointLocalRef(dm,kelem,ot,&rhs);
        VecGetArray(vec_rhs,&xx);
        for(unsigned kne=0; kne<NpE*_meqn; kne++)
            rhs[kne] = xx[kne];
        VecRestoreArray(vec_rhs,&xx);

    }
    DMRestoreLocalVector(dm, &locU);
    VecRestoreArray(out, &ot);

    isInfinityOrNAN(out, "NAN/INF encountered in RHS Vector of DG step-function");

    REAL newDt =  2./3.*_cfl*_geom->dtscale2D()*(_geom->rMin()/maxSpeed);
    status.setStatus(true);
    status.setSuggestedDt(newDt);

    // deallocate
    VecDestroy(&q_vol);
    VecDestroy(&q_surf);
    VecDestroy(&xcoord);
    VecDestroy(&ycoord);

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
