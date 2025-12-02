#include "aphwcalculategradients.h"
#include <wxmath.h>
#include "petsc_compat.h"  // PETSc API compatibility for version 3.19+

// Calculate the gradients using at the nodes for each element.
// This calculation is done in accordance with Chapter 7 of the
// Hesthaven and Warburton nodal DG book pg. 317.

template <typename REAL>
void
ApHWCalculateGradients<REAL>::setup(const WxCryptSet& wxc, DM dm)
{
  // call base class setup
  WxNodalDGLimiter<REAL>::setup(wxc, dm);

  _dm = dm;

  PetscInt kStart, kEnd;
  DMPlexGetHeightStratum(_dm, 0, &kStart, &kEnd);
  _Klocal = kEnd-kStart;

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

  _meqn = _eqnSet.totalEqns();

  // read list of BC subsolvers
  std::vector<WxAny> bcs;//, lbs;
  bcs = wxc.template get<std::vector<WxAny> >("boundaryConditions");

  std::vector<WxAny>::const_iterator i;
  for (i=bcs.begin(); i!=bcs.end(); ++i)
      _bcSubSolvers.push_back( wx_any_cast<std::string>(*i) );

  // Allocate memory for calculating the Elements and Face Gradients
  dVdxE1 = alloc_1d<REAL>(_meqn*_Klocal);
  dVdyE1 = alloc_1d<REAL>(_meqn*_Klocal);
  dVdxE2 = alloc_1d<REAL>(_meqn*_Klocal);
  dVdyE2 = alloc_1d<REAL>(_meqn*_Klocal);
  dVdxE3 = alloc_1d<REAL>(_meqn*_Klocal);
  dVdyE3 = alloc_1d<REAL>(_meqn*_Klocal);
  dVdxC0 = alloc_1d<REAL>(_meqn*_Klocal);
  dVdyC0 = alloc_1d<REAL>(_meqn*_Klocal);

}

template <typename REAL>
ApHWCalculateGradients<REAL>::~ApHWCalculateGradients()
{
    delete [] dVdxE1;
    delete [] dVdyE1;
    delete [] dVdxE2;
    delete [] dVdyE2;
    delete [] dVdxE3;
    delete [] dVdyE3;
    delete [] dVdxC0;
    delete [] dVdyC0;
    DMDestroy(&_dm);
}

template <typename REAL>
void
ApHWCalculateGradients<REAL>::calculateGradients(wxNodalDGgeometry2D<REAL> *geom, WxCubature2d<REAL> *cub, Vec qk, Vec q_grads)
{
    Vec local_out, local_in;
    // create local vector
    DMGetLocalVector(_dm, &local_in);
    DMGetLocalVector(_dm, &local_out);

    int kNodes = geom->NpElem();         // Number of nodes per elements
    int NpF = geom->NpFaces();           // Number of nodes per face

    // used to access data on the vectors
    const PetscScalar *u;
    PetscScalar *v, *qIn, *qOut;

    // get local values of the global vector in into locX
    DMGlobalToLocalBegin(_dm, qk, INSERT_VALUES, local_in);
    DMGlobalToLocalEnd(_dm, qk, INSERT_VALUES, local_in);
    VecGetArrayRead(local_in, &u);

    int connect[6]; // connectivity information element-to-element-to-edge
    REAL geoFacts[5], normals[9];

    // get id and values of nodes at faces
    int faceIDs[3*NpF];
    REAL QM[3*NpF*_meqn], QP[3*NpF*_meqn];
    geom->returnFmask(faceIDs);

    PetscInt kStart, kEnd, kEndInterior;
    DMPlexGetHeightStratum(_dm, 0, &kStart, &kEnd);
    DMPlexGetHybridBounds(_dm, &kEndInterior, NULL, NULL, NULL);

    for(int eNum=kStart; eNum<kEndInterior; eNum++)
    {
        // Get node values for this element
        DMPlexPointLocalRead(_dm, eNum, u, &qVal);
        for(unsigned kk=0; kk<kNodes*_meqn; kk++)
            q_vol[kk] = qVal[kk];

        // calculate edge normals and element Jacobian
        REAL geoFacts[5], normals[3*NfE];
        _geom->GeometricFactors2d(eNum,geoFacts); // [drdx, dsdx, drdy, dsdy, J]
        _geom->Normals2d(eNum,normals); // [nx_edge1,ny_edge1,length_edge1, nx_edge2, ny_edge2 ...]

        /** *******************************************************
         *  *******************************************************
         *  Evaluate Volume Integral and Sources
         *  *******************************************************
         *  *******************************************************
         */

        // interpolate nodes values into cubature points
        cub->interpolatedTOCubatures(_meqn,q_vol,Iq_vol);
        cub->interpolatedTOCubatures(1,xcoord,IXcoords);
        cub->interpolatedTOCubatures(1,ycoord,IYcoords);


    }


    VecRestoreArrayRead(local_in, &u);
    VecRestoreArray(local_out, &v);

    DMLocalToGlobalBegin(_dm, local_out, INSERT_VALUES, q_grads);
    DMLocalToGlobalEnd(_dm, local_out, INSERT_VALUES, q_grads);

    DMRestoreLocalVector(_dm, &local_in);
    DMRestoreLocalVector(_dm, &local_out);

    VecDestroy(&local_out);
    VecDestroy(&local_in);
}


template<typename REAL>
void
ApHWCalculateGradients<REAL>::applyBc(int bcNum, REAL *xc, REAL *nx, REAL *q, REAL *qaux, REAL *AreaInts, REAL *qBC)
{
    // apply boundary conditions
    ApSubSolver<REAL>* ss = this->getParent()->getSubSolver( _bcSubSolvers.at(bcNum-1) );
    // cast this to the a grid BC and call step function
    dynamic_cast<WxGridBC<REAL>* >(ss)->applyToArray(xc,nx,q,qaux,AreaInts,qBC);
}

template <typename REAL>
PetscErrorCode
ApHWCalculateGradients<REAL>::isInfinityOrNAN(Vec f, std::string location)
{
    PetscReal fnorm;
    VecNormBegin(f,NORM_2,&fnorm);	/* fnorm <- ||F||  */
    VecNormEnd(f,NORM_2,&fnorm);
    if (PetscIsInfOrNanReal(fnorm))
    {
        WxLogger *l = WxLogger::get("apollo-root.console");
        WxLogStream errStrm = l->getErrorStream();
        errStrm << location ;
        exit(1); // abort execution
    }
    return 0;
}

// instantiations
template class ApHWCalculateGradients<float>;
template class ApHWCalculateGradients<double>;
