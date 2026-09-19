#include "apdiscontinuitymitigation.h"
#include <wxmath.h>
#include "petsc_compat.h"  // PETSc API compatibility for version 3.19+

template <typename REAL>
void
ApDiscontinuityMitigation<REAL>::setup(const WxCryptSet& wxc, DM dm)
{
  // call base class setup
  WxNodalDGLimiter<REAL>::setup(wxc, dm);

  _dm = dm;

   _polyOrder = wxc.template get<int>("polynomialOrder");

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

  PetscInt kStart, kEnd;
  DMPlexGetHeightStratum(_dm, 0, &kStart, &kEnd);
  _Klocal = kEnd-kStart;

  // read list of Boundary Conditions
//  std::vector<WxAny> bcs;
//  bcs = wxc.template get<std::vector<WxAny> >("boundaryConditions");
//  std::vector<WxAny>::const_iterator i;
//  for (i=bcs.begin(); i!=bcs.end(); ++i)
//      _bcSubSolvers.push_back( wx_any_cast<std::string>(*i) );

  // read the list of equation system needed
//  std::vector<WxAny> eqns;
//  eqns = wxc.template get<std::vector<WxAny> >("equations");
//  for (i=eqns.begin(); i!=eqns.end(); ++i)
//      _bcSubSolvers.push_back( wx_any_cast<std::string>(*i) );

  _qM = alloc_1d<REAL>(_meqn);
  _qP = alloc_1d<REAL>(_meqn);
  _surfacesIntegral = alloc_1d<REAL>(_meqn);

}

template <typename REAL>
ApDiscontinuityMitigation<REAL>::~ApDiscontinuityMitigation()
{
    DMDestroy(&_dm);
    delete [] _qM;
    delete [] _qP;
    delete [] _surfacesIntegral;
}

template <typename REAL>
void
ApDiscontinuityMitigation<REAL>::detectElements(wxNodalDGgeometry2D<REAL> *geom, WxCubature2d<REAL> *cub, Vec qk, int *limElems)
{
    int kNodes = geom->NpElem();         // Number of nodes per element
    int NpF = geom->NpFaces();           // Number of nodes per face
    int NfE = geom->NfElem();            // Number of faces per element
    int Ncubature= cub->numCubaturePoints(); // Number of cubature points
    int Ngauss = cub->numGaussianPoints(); // Number of Gaussian points per edge/face

    // used to access data on the vectors
    PetscScalar *u, *qIn;
    PetscScalar *v, *qOut;

    PetscInt kStart, kEnd, kEndInterior;
    DMPlexGetHeightStratum(_dm, 0, &kStart, &kEnd);
    DMPlexGetHybridBounds(_dm, &kEndInterior, NULL, NULL, NULL);

    VecGetArray(qk, &u);

    // Coordinates
    REAL xcoord[kNodes], ycoord[kNodes];
    REAL xc[4]; xc[0]= this->getCurrentTime(); // coordinates
    REAL nx[2]; // normals
    int connect[2*NfE]; // connectivity information element-to-element-to-edge
    REAL Xcrd[NfE*Ngauss], Ycrd[NfE*Ngauss];;

    // Volume values
    REAL q_vol[kNodes*_meqn];

    // Surface intr
    REAL qtemp[kNodes*_meqn], QP[NfE*Ngauss*_meqn], QM[NfE*Ngauss*_meqn];
    REAL qgtemp[NfE*Ngauss*_meqn];

    for(int eNum=kStart; eNum<kEndInterior; eNum++)
    {
        // get coordinates of all nodes
        for(unsigned nodes=0; nodes<NpE; nodes++)
        {
            xcoord[nodes] = _geom->Xcoordinate(eNum,nodes);
            ycoord[nodes] = _geom->Ycoordinate(eNum,nodes);
        }

        // Get node values for this element
        DMPlexPointLocalRef(_dm, eNum, u, &qIn);
        for(unsigned kk=0; kk<kNodes*_meqn; kk++)
            q_vol[kk] = qIn[kk];

        // calculate edge normals and element Jacobian
        REAL normals[3*NfE];
        geom->Normals2d(eNum,normals); // [nx_edge1,ny_edge1,length_edge1, nx_edge2, ny_edge2 ...]

        // interpolate nodal values to surface Gassian quadrature points
        _cub->nodesTOSurfaceGaussians(_meqn,q_vol,QM);

        // Interpolate the coordinates at the element edges
        cub->nodesTOSurfaceGaussians(1,xcoord,Xcrd);
        cub->nodesTOSurfaceGaussians(1,ycoord,Ycrd);

        // Flux Gather
        // get the values at the Gaussian points of the adjacent elements
        geom->ElementTOElementANDFace(kelem,connect);
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

                    applyBc(abs(edgeNum), xc, nx, _qM, _qauxM, 0, _qP); // NOTE: boundary conditions that depend on the surface integral calculation will not be done properly on this step

                    for(unsigned comp=0; comp<_meqn; comp++)
                        QP[(edge*Ngauss+gpoint)*_meqn+comp] = _qP[comp];
                }
            }
            else
            {
                int plusElem = connect[2*edge];
                // get the values on the element adjacent to this edge
                DMPlexPointLocalRef(_dm, plusElem, u, &qIn);
                for(unsigned kk=0; kk<kNodes*_meqn; kk++)
                    qtemp[kk] = qIn[kk];

                // interpolate values at all edges of opposing element
                _cub->nodesTOSurfaceGaussians(_meqn,qtemp,qgtemp);

                // only use the Gaussian values of the needed edge
                for(unsigned gpoint=0; gpoint<Ngauss; gpoint++)
                    for(unsigned comp=0; comp<_meqn; comp++)
                        QP[(edge*Ngauss+gpoint)*_meqn+comp] = qgtemp[(edgeNum*Ngauss+Ngauss-1-gpoint)*_meqn+comp];
            }
        }

        // Substract values on both sides of the face
        for(unsigned subs=0; subs>NfE*Ngauss*_meqn; subs++)
            qgtemp[subs] = fabs(QM[subs] - QP[subs]);

        // Specify which variable to use as trigger
        for(unsigned kk=0; kk<Ngauss*NfE; kk++)
        {
            int curEdge = kk/Ngauss;

            for(unsigned comp=0; comp<_meqn; comp++)
                _qM[comp] = qgtemp[kk*_meqn+comp];

            _eqnSet.DGLimiterTrigger(_qM,_qP);

            for(unsigned comp=0; comp<_meqn; comp++){
                if(_qP[comp]!=0.0)
                    qgtemp[kk*_meqn+comp] = _qP[comp]/normals[3*curEdge+2]/QM[kk*_meqn+comp]/pow(geom->rMin(),(_polyOrder+1)/2.);
            }
        }

        // Discontinuity detector
        REAL I_j[_meqn];
        cub->discontinuityDetectorIntegral(qgtemp,I_j);
    }

    VecRestoreArray(qk, &u);

}

template <typename REAL>
void
ApDiscontinuityMitigation<REAL>::applyLimiter(wxNodalDGgeometry2D<REAL> *geom, WxCubature2d<REAL> *cub, Vec qk, Vec q_limited)
{
    int kNodes = geom->NpElem();         // Number of nodes per elements
    int NpF = geom->NpFaces();           // Number of nodes per face

}

template<typename REAL>
void
ApDiscontinuityMitigation<REAL>::applyBc(int bcNum, REAL *xc, REAL *nx, REAL *q, REAL *qaux, REAL *AreaInts, REAL *qBC)
{
    // apply boundary conditions
    ApSubSolver<REAL>* ss = this->getParent()->getSubSolver( _bcSubSolvers.at(bcNum-1) );
    // cast this to the a grid BC and call step function
    dynamic_cast<WxGridBC<REAL>* >(ss)->applyToArray(xc,nx,q,qaux,AreaInts,qBC);
}

// instantiations
//template class ApDiscontinuityMitigation<float>;
template class ApDiscontinuityMitigation<double>;
