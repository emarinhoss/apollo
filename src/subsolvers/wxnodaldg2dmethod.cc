#include "wxnodaldg2dmethod.h"

// WarpX lib includes
#include <wxcreator.h>
#include <wxlogger.h>
#include <wxlogstream.h>

# include <wxmpimsg.h>

#ifdef _OPENMP
#include <omp.h>
#endif

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
    delete [] _areaInts;
    delete [] _AgregateAreaIntegral;
    delete [] _df;
    delete [] _amdq;
    delete [] _apdq;
    delete [] _sx;
    delete [] _sy;
    delete [] _qauxM;
    delete [] _qauxP;
    delete [] _numericalFLux;
    //delete [] _filterdiag;
    //delete [] _filterMatrix;
    free_2d_c(_wave,_meqn,_mwave);
    delete _geom;
    delete _cub;
    delete _initFunc;
//    VecDestroy(&locU);
//    DMDestroy(&_dm);
}

// M - denotes the interior of the edge
// P - denotes the exterior of the edge

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

  _areaSet.setNumEqns(_meqn); // set no of equations
  _areaSet.setup(wxc);

  _qM = alloc_1d<REAL>(_meqn);
  _qP = alloc_1d<REAL>(_meqn);
  _qauxM = alloc_1d<REAL>(_meqn);
  _qauxP = alloc_1d<REAL>(_meqn);
  _df = alloc_1d<REAL>(_meqn); // jump
  _src = alloc_1d<REAL>(_meqn);
  _areaInts = alloc_1d<REAL>(_meqn);
  _AgregateAreaIntegral = alloc_1d<REAL>(_meqn);
  _fM = alloc_1d<REAL>(_meqn);
  _fP = alloc_1d<REAL>(_meqn);
  _gM = alloc_1d<REAL>(_meqn);
  _gP = alloc_1d<REAL>(_meqn);
  _numericalFLux = alloc_1d<REAL>(_meqn);
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

  // read list of BC subsolvers
  std::vector<WxAny> bcs;//, lbs;
  bcs = wxc.template get<std::vector<WxAny> >("boundaryConditions");

  std::vector<WxAny>::const_iterator i;
  for (i=bcs.begin(); i!=bcs.end(); ++i)
      _bcSubSolvers.push_back( wx_any_cast<std::string>(*i) );

  _haveLimiter = false;
  if (wxc.has("Limiter"))
  {
      std::vector<WxAny> lmt;
      lmt = wxc.template get<std::vector<WxAny> >("Limiter");
      _haveLimiter = true;
      for (i=lmt.begin(); i!=lmt.end(); ++i)
          _limiterSubSolvers.push_back( wx_any_cast<std::string>(*i) );
  }

  _calculateGradients = false;
  if (wxc.has("Gradients"))
  {
      std::vector<WxAny> grt;
      grt = wxc.template get<std::vector<WxAny> >("Gradients");
      _calculateGradients = true;
      for (i=grt.begin(); i!=grt.end(); ++i)
          _gradientSubSolvers.push_back( wx_any_cast<std::string>(*i) );
  }

  // Calculate connectivity, coordinates and Matrices
  _geom = new wxNodalDGgeometry2D<REAL>(_dm, _meqn, _polyOrder);

  // Evaluate cubatures, surface Gaussian points and derivatives
  _cub = new WxCubature2d<REAL>(_dm,_meqn,_polyOrder,_geom->inverseVandermonde());

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

    for (k = kStart; k < kEndInterior; ++k)
    {
        for(unsigned node=0; node<_geom->NpElem(); node++){
            txo[1] = _geom->Xcoordinate(k,node);
            txo[2] = _geom->Ycoordinate(k,node);

            _initFunc->func(3, txo, d);
            PetscScalar *xc;

            // reference this cell to the proper location on the solution
            // vector
            DMPlexPointGlobalRef(_dm,k,x,&xc);
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

    isInfinityOrNAN(out, "NAN/INF in initialization!\n");

    // Suggested initial dt
    REAL timeStep = 2./3.*_cfl*_geom->dtscale2D()*(_geom->rMin()/maxSpeed);

    REAL SmallestDT;
    PetscBarrier((PetscObject) _dm);
    MPI_Allreduce(&timeStep, &SmallestDT, 1,
                      MPI_DOUBLE, MPI_MIN,
                      MPI_COMM_WORLD);

    this->setDt(SmallestDT);

    for(unsigned kx=0; kx<_meqn; kx++)
        _AgregateAreaIntegral[kx] = 0.0;

}

template <typename REAL>
WxStepperStatus<REAL>
WxNodalDG2dMethod<REAL>::step(REAL t, REAL dt, Vec in, Vec out)
{
    DM dataManage;
    Vec local_out, local_in, gradients;

    // Apply Limiter
    if(_haveLimiter){
//      isInfinityOrNAN(in, "NAN/INF encountered before limiting occurs.\n");
        applyLimiter(in,in);
        isInfinityOrNAN(in, "NAN/INF encountered in limiter vector of DG step-function.\n");
    }

    if(_calculateGradients){
        VecDuplicate(in,&gradients);
        calculateGradients(in,gradients);
        isInfinityOrNAN(in, "NAN/INF encountered in gradient vector of DG step-function.\n");
    }

    WxStepperStatus<REAL> status;

    const PetscScalar *u;
    PetscScalar *ot, *rhs;
    REAL maxSpeed=0.0;

    // get data Manager
    VecGetDM(in,&dataManage);
    // create local vector
    DMGetLocalVector(dataManage, &local_in);
    DMGetLocalVector(dataManage, &local_out);

    // zero entries of the vectors that will be used to store
    // information
//    VecZeroEntries(local_in);
//    VecZeroEntries(local_out);

    // get local values of the global vector in into locX
    DMGlobalToLocalBegin(dataManage, in, INSERT_VALUES, local_in);
    DMGlobalToLocalEnd(dataManage, in, INSERT_VALUES, local_in);

//    isInfinityOrNAN(local_in, "NAN/INF in input Vector to DG step-function");

    // get first and last element number
    PetscInt kStart, kEnd, kEndInterior;
    DMPlexGetHeightStratum(dataManage, 0, &kStart, &kEnd);
    DMPlexGetHybridBounds(dataManage, &kEndInterior, NULL , NULL, NULL);
    VecGetArrayRead(local_in, &u);
    VecGetArray(local_out, &ot);

    int NpE = _geom->NpElem(); // Number of nodes per Element
    int NfE = _geom->NfElem(); // Number of edge/faces per Element
    int Ncubature= _cub->numCubaturePoints(); // Number of cubature points
    int Ngauss = _cub->numGaussianPoints(); // Number of Gaussian points per edge/face

    // Total area integrals (accumulated from all elements)
    REAL TotalAreaInt[_meqn];
    for(unsigned eqs=0; eqs<_meqn; eqs++)
        TotalAreaInt[eqs] = 0.0;

    // OpenMP parallelization of element loop
    // Each element computation is independent, making this embarrassingly parallel
    #pragma omp parallel for reduction(max:maxSpeed) schedule(static)
    for(unsigned kelem=kStart; kelem<kEndInterior; kelem++)
    {
        // Thread-local temporary arrays (all arrays are now private to each thread)
        REAL xcoord[NpE], ycoord[NpE];
        REAL xc[5]; xc[0]=t; xc[4] = dt;
        REAL nx[2];
        int connect[2*NfE];

        REAL q_vol[NpE*_meqn], vec_rhs[NpE*_meqn], volInt[NpE*_meqn];
        REAL Iq_vol[Ncubature*_meqn], If_vol[Ncubature*_meqn], Ig_vol[Ncubature*_meqn],
                ISrc_vol[Ncubature*_meqn], Iarea_vol[Ncubature*_meqn];
        REAL IXcoords[Ncubature], IYcoords[Ncubature];

        REAL qtemp[NpE*_meqn], q_surf[NpE*_meqn];
        REAL QP[NfE*Ngauss*_meqn], QM[NfE*Ngauss*_meqn], numFlux[NfE*Ngauss*_meqn],
                Xcrd[NfE*Ngauss], Ycrd[NfE*Ngauss];
        REAL qgtemp[NfE*Ngauss*_meqn];

        REAL Qvar[_meqn], Qvaraux[_meqn], Fflux[_meqn], Gflux[_meqn], AreaIntegrals[_meqn];
        REAL geoFacts[5], normals[3*NfE];

        PetscScalar *qVal, *rhs;

        // get coordinates of all nodes
        for(unsigned nodes=0; nodes<NpE; nodes++)
        {
            xcoord[nodes] = _geom->Xcoordinate(kelem,nodes);
            ycoord[nodes] = _geom->Ycoordinate(kelem,nodes);
        }

        // Get node values for this element
        DMPlexPointLocalRead(dataManage, kelem, u, &qVal);
        for(unsigned kk=0; kk<NpE*_meqn; kk++)
            q_vol[kk] = qVal[kk];

        // calculate edge normals and element Jacobian
        REAL geoFacts[5], normals[3*NfE];
        _geom->GeometricFactors2d(kelem,geoFacts); // [drdx, dsdx, drdy, dsdy, J]
        _geom->Normals2d(kelem,normals); // [nx_edge1,ny_edge1,length_edge1, nx_edge2, ny_edge2 ...]

        /** *******************************************************
         *  *******************************************************
         *  Evaluate Volume Integral and Sources
         *  *******************************************************
         *  *******************************************************
         */

        // interpolate nodes values into cubature points
        _cub->interpolatedTOCubatures(_meqn,q_vol,Iq_vol);
        _cub->interpolatedTOCubatures(1,xcoord,IXcoords);
        _cub->interpolatedTOCubatures(1,ycoord,IYcoords);

        for(unsigned point=0; point<Ncubature; point++)
        {
            // Coordinates
            xc[1] = IXcoords[point];
            xc[2] = IYcoords[point];

            // Evaluate Flux at each cubature point
            for(unsigned component=0; component<_meqn; component++)
                Qvar[component] = Iq_vol[point*_meqn+component];

//            REAL AA = Qvar[15];
            _eqnSet.flux(0, xc, Qvar, Qvaraux, Fflux);
            _eqnSet.flux(1, xc, Qvar, Qvaraux, Gflux);
            _srcSet.sourceTerms(xc, Qvar, Qvaraux, _src);
            _areaSet.areaTerms(xc, Qvar, Qvaraux, _areaInts);

            for(unsigned component=0; component<_meqn; component++)
            {
                  If_vol[point*_meqn+component] = Fflux[component];
                  Ig_vol[point*_meqn+component] = Gflux[component];
                ISrc_vol[point*_meqn+component] = _src[component];
               Iarea_vol[point*_meqn+component] = _areaInts[component];
            }
        }

        // Compute volume terms, which includes the sources: (dphidx, F) + (dphidy, G) + (phi, S)
        _cub->evaluatedVolumeIntegrals(xcoord,ycoord,If_vol,Ig_vol,ISrc_vol,volInt);

        // Calculate the area integrals
        _cub->CalculateAreaIntegrals(xcoord,ycoord,Iarea_vol,AreaIntegrals);

        /** *******************************************************
         *  *******************************************************
         *  Evaluate Surface Integral
         *  *******************************************************
         *  *******************************************************
         */

        // interpolate nodal values to surface Gassian quadrature points
        _cub->nodesTOSurfaceGaussians(_meqn,q_vol,QM);

        // Interpolate the coordinates at the element edges
        _cub->nodesTOSurfaceGaussians(1,xcoord,Xcrd);
        _cub->nodesTOSurfaceGaussians(1,ycoord,Ycrd);

        // Flux Gather
        // get the values at the Gaussian points of the adjacent elements
        _geom->ElementTOElementANDFace(kelem,connect);
        for(unsigned edge=0; edge<NfE; edge++)
        {
            // gather only the values for the needed edge
            int edgeNum = connect[2*edge+1];
            // If edgeNum is negative this edge is a physical boundary.
            if(edgeNum<0)
            {
                // Apply Boundary conditions
                for(unsigned gpoint=0; gpoint<Ngauss; gpoint++)
                {
                    xc[1] = Xcrd[edge*Ngauss+gpoint];
                    xc[2] = Ycrd[edge*Ngauss+gpoint];

                    nx[0] = normals[3*edge];
                    nx[1] = normals[3*edge+1];

                    for(unsigned cmp=0; cmp<_meqn; cmp++)
                        _qM[cmp] = QM[(edge*Ngauss+gpoint)*_meqn+cmp];

                    applyBc(abs(edgeNum), xc, nx, _qM, _qauxM, _AgregateAreaIntegral, _qP);

                    for(unsigned comp=0; comp<_meqn; comp++)
                        QP[(edge*Ngauss+gpoint)*_meqn+comp] = _qP[comp];
                }
            }
            else
            {
                int plusElem = connect[2*edge];
                // get the values on the element adjacent to this edge
                DMPlexPointLocalRead(dataManage, plusElem, u, &qVal);
                for(unsigned kk=0; kk<NpE*_meqn; kk++)
                    qtemp[kk] = qVal[kk];

                // interpolate values at all edges of opposing element
                _cub->nodesTOSurfaceGaussians(_meqn,qtemp,qgtemp);

                // only use the Gaussian values of the needed edge
                for(unsigned gpoint=0; gpoint<Ngauss; gpoint++)
                    for(unsigned comp=0; comp<_meqn; comp++)
                        QP[(edge*Ngauss+gpoint)*_meqn+comp] = qgtemp[(edgeNum*Ngauss+Ngauss-1-gpoint)*_meqn+comp];
            }
        }

        // Calculate the Numerical Flux at each
        // Gaussian quadrature point at each edge
        for(unsigned kk=0; kk<Ngauss*NfE; kk++)
        {
            int curEdge = kk/Ngauss;

            for(unsigned comp=0; comp<_meqn; comp++){
                _qM[comp] = QM[kk*_meqn+comp];
                _qP[comp] = QP[kk*_meqn+comp];}

            REAL lambda; // Fastest propagating wave speed

            // face normals for this edge
            nx[0] = normals[3*curEdge+0];
            nx[1] = normals[3*curEdge+1];

            // Evaluate the numerical flux and get the speed of the
            // fastest propagating wave
            _eqnSet.DGnumericalFlux(nx,_qM,_qP,_numericalFLux,&lambda);
            lambda = sqrt(lambda*lambda);

            // find maximum propagation speed in the entire domain;
            // this will be used to adjust the timestep
            maxSpeed = dmax(maxSpeed,lambda);

            for(unsigned comp=0; comp<_meqn; comp++)
                numFlux[kk*_meqn+comp] = _numericalFLux[comp]*normals[3*curEdge+2];
        }

        // Compute surface integral contribution to all nodes
        _cub->calculateSurfaceIntegral(numFlux,q_surf);

        // add surface and volume contributions
        for(unsigned kk=0; kk<NpE*_meqn; kk++){
            volInt[kk] -= q_surf[kk];

            if(q_surf[kk]!=q_surf[kk]){
                WxLogger::get("apollo-root.console")->
                  error("*** NaN in surface integral evaluation RHS ***\n");
                exit(1); // abort execution
            }

            if(volInt[kk]!=volInt[kk]){
                WxLogger::get("apollo-root.console")->
                  error("*** NaN in volume integral evaluation RHS ***\n");
                exit(1); // abort execution
            }
        }

        // Multiply by the inverse Mass Matrix
        _geom->multiplyBYinverseMassMatrix(volInt,vec_rhs);

        DMPlexPointLocalRef(dataManage,kelem,ot,&rhs);
        for(unsigned kne=0; kne<NpE*_meqn; kne++){
            rhs[kne] = vec_rhs[kne]/geoFacts[4];

            if(vec_rhs[kne]!=vec_rhs[kne]){
                WxLogger::get("apollo-root.console")->
                  error("*** NaN after Mass matrix multiplication RHS ***\n");
                exit(1); // abort execution
            }

        }

        // compute \int q\cdot dA, add contribution from all elements
//        REAL elementArea = _geom->elementArea(kelem);
        #pragma omp critical
        {
            for(unsigned ar=0; ar<_meqn; ar++)
                TotalAreaInt[ar] += AreaIntegrals[ar];
        }
    }

    VecRestoreArrayRead(local_in, &u);
    VecRestoreArray(local_out, &ot);

//    isInfinityOrNAN(local_out, "NAN/INF encountered in RHS Vector of DG step-function.\n");

    DMLocalToGlobalBegin(dataManage, local_out, INSERT_VALUES, out);
    DMLocalToGlobalEnd(dataManage, local_out, INSERT_VALUES, out);

    DMRestoreLocalVector(dataManage, &local_in);
    DMRestoreLocalVector(dataManage, &local_out);

    REAL newDt =  2./3.*_cfl*_geom->dtscale2D()*(_geom->rMin()/maxSpeed);
    status.setStatus(true);
    status.setSuggestedDt(newDt);

    // add the surface integral contributions from all processors
    REAL inValue, outValue;
//    REAL AB = _AgregateAreaIntegral[15];
    for(unsigned numeq=0; numeq<_meqn; numeq++)
    {
        inValue = TotalAreaInt[numeq];
        MPI_Allreduce(&inValue, &outValue, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        _AgregateAreaIntegral[numeq] = outValue;
    }

//    REAL AC = _AgregateAreaIntegral[15];
//    REAL AD = AC-AB;
    VecDestroy(&local_out);
    VecDestroy(&local_in);
    return status;
}

template <typename REAL>
PetscErrorCode
WxNodalDG2dMethod<REAL>::isInfinityOrNAN(Vec f, std::string location)
{
    PetscReal fnorm;
    VecNormBegin(f,NORM_2,&fnorm);	/* fnorm <- ||F||  */
    VecNormEnd(f,NORM_2,&fnorm);
    if (PetscIsInfOrNanReal(fnorm))
    {
        //VecView(f,PETSC_VIEWER_STDOUT_WORLD);
        //REAL test = 0.0;
        WxLogger *l = WxLogger::get("apollo-root.console");
        WxLogStream errStrm = l->getErrorStream();
        errStrm << location ;
//        PetscFinalize();
        exit(1); // abort execution
    }
    return 0;
}

template<typename REAL>
void
WxNodalDG2dMethod<REAL>::applyBc(int bcNum, REAL *xc, REAL *nx, REAL *q, REAL *qaux, REAL *AreaInts, REAL *qBC)
{
    // apply boundary conditions
    ApSubSolver<REAL>* ss = this->getParent()->getSubSolver( _bcSubSolvers.at(bcNum-1) );
    // cast this to the a grid BC and call step function
    dynamic_cast<WxGridBC<REAL>* >(ss)->applyToArray(xc,nx,q,qaux,AreaInts,qBC);
}

template<typename REAL>
void
WxNodalDG2dMethod<REAL>::applyLimiter(Vec Qin, Vec Qlimited)
{
    // apply limiters
    ApSubSolver<REAL>* ss = this->getParent()->getSubSolver( _limiterSubSolvers.at(0));
    // cast this to the limiter and call step function
    dynamic_cast<WxNodalDGLimiter<REAL>* >(ss)->applyToVector(_geom,_cub,Qin,Qlimited);
}

template<typename REAL>
void
WxNodalDG2dMethod<REAL>::calculateGradients(Vec Qin, Vec Qgrads)
{
    // apply limiters
    ApSubSolver<REAL>* ss = this->getParent()->getSubSolver( _gradientSubSolvers.at(0));
    // cast this to the limiter and call step function
    dynamic_cast<ApNodalDGcalculateGradients<REAL>* >(ss)->calculateGradients(_geom,_cub,Qin,Qgrads);
}

// instantiations
template class WxNodalDG2dMethod<float>;
template class WxNodalDG2dMethod<double>;
