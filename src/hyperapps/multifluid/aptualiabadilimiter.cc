#include "aptualiabadilimiter.h"
#include <wxmath.h>
#include "petsc_compat.h"  // PETSc API compatibility for version 3.19+

// Limit the Euler solution using slope limiting adapted from
// A SLOPE LIMITING PROCEDURE IN DISCONTINUOUS GALERKIN FINITE ELEMENT METHOD FOR
// GASDYNAMICS APPLICATIONS. SHUANGZHANG TU AND SHAHROUZ ALIABADI
// INTERNATIONAL JOURNAL OF NUMERICAL ANALYSIS AND MODELING, Volume 2, Number 2, Pages 163

template <typename REAL>
void
WxTuAliabadiLimiter<REAL>::setup(const WxCryptSet& wxc, DM dm)
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
WxTuAliabadiLimiter<REAL>::~WxTuAliabadiLimiter()
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
WxTuAliabadiLimiter<REAL>::applyLimiter(wxNodalDGgeometry2D<REAL> *geom, WxCubature2d<REAL> *cub, Vec qk, Vec q_limited)
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

    // get average vectors and matrices
    REAL AVE[kNodes], dropAVE[kNodes*kNodes];
    geom->LimiterElementAVE(AVE,dropAVE);
    REAL lim_dx[kNodes], lim_dy[kNodes];

    int connect[6]; // connectivity information element-to-element-to-edge
    REAL geoFacts[5], normals[9];

    // get id and values of nodes at faces
    int faceIDs[3*NpF];
    REAL QM[3*NpF*_meqn], QP[3*NpF*_meqn];
    geom->returnFmask(faceIDs);

    // array that has the cell conservative/primitive averages for every element
    REAL cellAVEcons[_meqn*_Klocal], cellAVEprim[_meqn*_Klocal];

    PetscInt kStart, kEnd, kEndInterior;
    DMPlexGetHeightStratum(_dm, 0, &kStart, &kEnd);
    DMPlexGetHybridBounds(_dm, &kEndInterior, NULL, NULL, NULL);

    for(int eNum=kStart; eNum<kEndInterior; eNum++)
    {
        // Not every cell in the stratum is a triangle; see
        // wxNodalDGgeometry2D::isTriangle. Asking a degenerate one for its
        // normals gives a zero-length edge and stops the run.
        if (!geom->isTriangle(eNum)) continue;

        // Cell centers
        REAL xc[4]={0,0,0,0}, yc[4]={0,0,0,0};
        // weights for face gradients
        REAL A0[4]={0,0,0,0};

        /**
          * Step 1: Compute geometric information for 4 element patch
          */

        // get neighbors ids
        geom->ElementTOElementANDFace(eNum,connect);
        // calculate edge normals
        geom->Normals2d(eNum,normals); // returns [nx_edge1,ny_edge1,length_edge1, nx_edge2, ny_edge2 ...]

        // calculate the element center
        geom->GeometricFactors2d(eNum,geoFacts); // returns [drdx, dsdx, drdy, dsdy, J]

        for(unsigned kk=0; kk<3; kk++)
        {
            xc[0] += AVE[kk]*geom->Xcoordinate(eNum,faceIDs[kk*NpF]);
            yc[0] += AVE[kk]*geom->Ycoordinate(eNum,faceIDs[kk*NpF]);
            A0[0] += 2./3.*AVE[kk]*geoFacts[4];
        }

        for(unsigned face=0; face<3; face++)
        {
            if(connect[2*face]<0){
                // for boundary edges mirror the cell center
                REAL H1 = 2.0*A0[0]/normals[3*face+2];
                xc[face+1] = xc[0] + 2.0*normals[3*face]*H1;
                yc[face+1] = yc[0] + 2.0*normals[3*face+1]*H1;
            }
            else
            {
                geom->GeometricFactors2d(connect[2*face],geoFacts);
                for(unsigned node=0; node<3; node++){
                    xc[face+1] += AVE[node]*geom->Xcoordinate(connect[2*face],faceIDs[node*NpF]);
                    yc[face+1] += AVE[node]*geom->Ycoordinate(connect[2*face],faceIDs[node*NpF]);
                    A0[face+1] += 2./3.*AVE[node]*geoFacts[4];}
            }
        }

        for(unsigned face=0; face<3; face++)
            A0[face+1] += A0[0];

        /**
          * Step 2: Find cell averages of conserved &
          * primitive variables in each 4 element patch
          */
        // extract fields from Q
        REAL ConsAve[_meqn][4], PrimAve[_meqn][4];
        REAL qCons[_meqn], qPrim[_meqn], qBC[_meqn];

        // for patch element 0
        DMPlexPointLocalRead(_dm, eNum, u, &qIn);
        _eqnSet.computeConservedAndPrimitiveAVEVariables(kNodes,qIn,AVE,qCons,qPrim);

        for(unsigned kk=0; kk<_meqn; kk++){
            ConsAve[kk][0] = qCons[kk];
            PrimAve[kk][0] = qPrim[kk];
            cellAVEprim[eNum*_meqn+kk] = qPrim[kk];
            cellAVEcons[eNum*_meqn+kk] = qCons[kk];
        }

        //
        for(unsigned edge=0; edge<3; edge++)
            for(unsigned nodes=0; nodes<NpF; nodes++)
                for(unsigned comp=0; comp<_meqn; comp++)
                    QM[(edge*NpF+nodes)*_meqn+comp] = qIn[faceIDs[edge*NpF+nodes]*_meqn+comp];

        // for patch elements 1-3
        for(unsigned elem=0; elem<3; elem++)
        {
            if(connect[2*elem]<0)
            {
                // face normal
                REAL NX[2], XC[2];

                // face normal
                NX[0] = normals[3*elem]; NX[1] = normals[3*elem+1];

                // cell center coordinates
                XC[0] = xc[elem+1]; XC[1] = yc[elem+1];

                for(unsigned kk=0; kk<_meqn; kk++)
                    qCons[kk] = ConsAve[kk][0];

                // Apply BC to cell averages of ghost cells
                applyBc(abs(connect[2*elem+1]),XC,NX,qCons,0,0,qBC);
                // Calculate Primitives
                _eqnSet.getPrimitiveVariable(qBC,qPrim);

                for(unsigned kk=0; kk<_meqn; kk++){
                    ConsAve[kk][elem+1] = qBC[kk];
                    PrimAve[kk][elem+1] = qPrim[kk];
                }

                // Apply BC to face nodes
                for(unsigned nodes=0; nodes<NpF; nodes++)
                {
                    // node coordinates
                    XC[0] = geom->Xcoordinate(eNum,faceIDs[elem*NpF+nodes]);
                    XC[1] = geom->Ycoordinate(eNum,faceIDs[elem*NpF+nodes]);

                    for(unsigned kk=0; kk<_meqn; kk++)
                        qCons[kk] = qIn[faceIDs[elem*NpF+nodes]*_meqn+kk];

                    // Apply BC
                    applyBc(abs(connect[2*elem+1]),XC,NX,qCons,0,0,qBC);

                    // store BC evaluation
                    for(unsigned kk=0; kk<_meqn; kk++)
                        QP[(elem*NpF+nodes)*_meqn+kk] = qBC[kk];
                }
            }
            else
            {
                DMPlexPointLocalRead(_dm, connect[2*elem], u, &qIn);
                _eqnSet.computeConservedAndPrimitiveAVEVariables(kNodes,qIn,AVE,qCons,qPrim);

                for(unsigned kk=0; kk<_meqn; kk++){
                    ConsAve[kk][elem+1] = qCons[kk];
                    PrimAve[kk][elem+1] = qPrim[kk];}

                // get adjacent node values
                int faceNum = connect[2*elem+1];
                for(unsigned nodes=0; nodes<NpF; nodes++)
                    for(unsigned kk=0; kk<_meqn; kk++)
                        QP[(elem*NpF+nodes)*_meqn+kk] = qIn[faceIDs[faceNum*NpF+NpF-1-nodes]*_meqn+kk];
            }
        }

        /**
         * Step 3: Compute average of primitive variables at face nodes
         */
        // Conserved variables face averages
        REAL fConsA[3*NpF*_meqn];

        for(unsigned kk=0; kk<3*NpF*_meqn; kk++)
            fConsA[kk] = 0.5*(QP[kk]+QM[kk]);

        // Get the primitive values now.
        // Only the values of the nodes at both ends of a
        // edge/face are needed.
        int endIDs[6]={faceIDs[0],faceIDs[NpF-1],faceIDs[NpF],faceIDs[2*NpF-1],faceIDs[2*NpF],faceIDs[3*NpF-1]};
        REAL fPrimA[6*_meqn], Conserved[_meqn], Primitive[_meqn];

        for(unsigned nodes=0; nodes<6; nodes++)
        {
            for(unsigned component=0; component<_meqn; component++)
                Conserved[component] = fConsA[endIDs[nodes]*_meqn+component];

            _eqnSet.getPrimitiveVariable(Conserved,Primitive);

            for(unsigned component=0; component<_meqn; component++)
                fPrimA[nodes*_meqn+component] = Primitive[component];

        }

        /** Step 4: Apply limiting procedure to each of the primitive variables */

        // Compute face gradients
        REAL xv1 = geom->Xcoordinate(eNum,endIDs[0]); REAL yv1 = geom->Ycoordinate(eNum,endIDs[0]);
        REAL xv2 = geom->Xcoordinate(eNum,endIDs[2]); REAL yv2 = geom->Ycoordinate(eNum,endIDs[2]);
        REAL xv3 = geom->Xcoordinate(eNum,endIDs[4]); REAL yv3 = geom->Ycoordinate(eNum,endIDs[4]);
        for(unsigned comp=0; comp<_meqn; comp++)
        {
            dVdxE1[_meqn*eNum+comp] =  0.5*((PrimAve[comp][1]-PrimAve[comp][0])*(yv2-yv1)
                                          + (fPrimA[0*_meqn+comp]-fPrimA[1*_meqn+comp])*(yc[1]-yc[0]) )/A0[1];

            dVdyE1[_meqn*eNum+comp] = -0.5*((PrimAve[comp][1]-PrimAve[comp][0])*(xv2-xv1)
                                          + (fPrimA[0*_meqn+comp]-fPrimA[1*_meqn+comp])*(xc[1]-xc[0]) )/A0[1];

            dVdxE2[_meqn*eNum+comp] =  0.5*((PrimAve[comp][2]-PrimAve[comp][0])*(yv3-yv2)
                                          + (fPrimA[2*_meqn+comp]-fPrimA[3*_meqn+comp])*(yc[2]-yc[0]) )/A0[2];

            dVdyE2[_meqn*eNum+comp] = -0.5*((PrimAve[comp][2]-PrimAve[comp][0])*(xv3-xv2)
                                          + (fPrimA[2*_meqn+comp]-fPrimA[3*_meqn+comp])*(xc[2]-xc[0]) )/A0[2];

            dVdxE3[_meqn*eNum+comp] =  0.5*((PrimAve[comp][3]-PrimAve[comp][0])*(yv3-yv2)
                                          + (fPrimA[4*_meqn+comp]-fPrimA[5*_meqn+comp])*(yc[3]-yc[0]) )/A0[3];

            dVdyE3[_meqn*eNum+comp] = -0.5*((PrimAve[comp][3]-PrimAve[comp][0])*(xv3-xv2)
                                          + (fPrimA[4*_meqn+comp]-fPrimA[5*_meqn+comp])*(xc[3]-xc[0]) )/A0[3];

            dVdxC0[_meqn*eNum+comp] = (A0[1]*dVdxE1[_meqn*eNum+comp] + A0[2]*dVdxE2[_meqn*eNum+comp]
                    + A0[3]*dVdxE3[_meqn*eNum+comp]) / (A0[1]+A0[2]+A0[3]);

            dVdyC0[_meqn*eNum+comp] = (A0[1]*dVdyE1[_meqn*eNum+comp] + A0[2]*dVdyE2[_meqn*eNum+comp]
                    + A0[3]*dVdyE3[_meqn*eNum+comp]) / (A0[1]+A0[2]+A0[3]);
        }
    }
    VecRestoreArrayRead(local_in, &u);

    // Cell center gradients of adjacent cells
    REAL dVdxC[3], dVdyC[3];

    // limited gradient nodes of each element
    REAL qdv[_meqn*kNodes];

    VecGetArray(local_out, &v);

    for(unsigned eNum=kStart; eNum<kEndInterior; eNum++)
    {
//        bool update = false;
        if (!geom->isTriangle(eNum)) continue;

        // get neighbors ids
        geom->ElementTOElementANDFace(eNum,connect);

        for(unsigned i=0; i<kNodes; i++){
            lim_dx[i] = 0.0; lim_dy[i] = 0.0;}

        for(unsigned i=0; i<kNodes; i++)
            for(unsigned j=0; j<kNodes; j++)
            {
                lim_dx[i] += dropAVE[i*kNodes+j]*geom->Xcoordinate(eNum,j);
                lim_dy[i] += dropAVE[i*kNodes+j]*geom->Ycoordinate(eNum,j);
            }

        // Loop over the components
        for(unsigned comp=0; comp<_meqn; comp++)
        {
            for(unsigned edge=0; edge<3; edge++)
            {
                // if this edge/face in the boundary, use egde/face gradients instead of cell center gradients
                if(connect[2*edge]<0)
                {
                    switch (edge) {
                    case 0:
                        dVdxC[edge] = dVdxE1[eNum*_meqn+comp];
                        dVdyC[edge] = dVdyE1[eNum*_meqn+comp];
                        break;
                    case 1:
                        dVdxC[edge] = dVdxE2[eNum*_meqn+comp];
                        dVdyC[edge] = dVdyE2[eNum*_meqn+comp];
                        break;
                    case 2:
                        dVdxC[edge] = dVdxE3[eNum*_meqn+comp];
                        dVdyC[edge] = dVdyE3[eNum*_meqn+comp];
                        break;
                    }
                }
                else
                {
                    int adjElement = connect[2*edge];
                    dVdxC[edge] = dVdxC0[adjElement*_meqn+comp];
                    dVdyC[edge] = dVdyC0[adjElement*_meqn+comp];
                }
            }

            // Build weights used in limiting
            REAL g1 = (dVdxC[0]*dVdxC[0] + dVdyC[0]*dVdyC[0]);
            REAL g2 = (dVdxC[1]*dVdxC[1] + dVdyC[1]*dVdyC[1]);
            REAL g3 = (dVdxC[2]*dVdxC[2] + dVdyC[2]*dVdyC[2]);

            REAL fac   = g1*g1+g2*g2+g3*g3;
            REAL epse  = 1.e-14;
            REAL fac3e = fac+3.0*epse;

            REAL w1 = (g2*g3 + epse) / fac3e;
            REAL w2 = (g1*g3 + epse) / fac3e;
            REAL w3 = (g1*g2 + epse) / fac3e;

            // Limit gradients
            REAL LdVdxC0 = w1*dVdxC[0] + w2*dVdxC[1] + w3*dVdxC[2];
            REAL LdVdyC0 = w1*dVdyC[0] + w2*dVdyC[1] + w3*dVdyC[2];

            for(unsigned nodes=0; nodes<kNodes; nodes++)
                qdv[nodes*_meqn+comp] = lim_dx[nodes]*LdVdxC0 + lim_dy[nodes]*LdVdyC0;
        }

        /** Step 5: Reconstruct conserved variables using cell averages and limited gradients  */

        DMPlexPointLocalRef(_dm, eNum, v, &qOut);

        for(unsigned nodes=0; nodes<kNodes; nodes++)
        {
            REAL avgPrimValues[_meqn], avgConsValues[_meqn], dValues[_meqn], limititedValues[_meqn];

            for(unsigned components=0; components<_meqn; components++)
            {
                avgPrimValues[components] = cellAVEprim[eNum*_meqn+components];
                avgConsValues[components] = cellAVEcons[eNum*_meqn+components];
                dValues[components] = qdv[nodes*_meqn+components];
            }

            // call the limited values from the eqn_system
            _eqnSet.tuAndAliabadiLimiter(avgConsValues,avgPrimValues,dValues,limititedValues);

            // limited values
            for(unsigned components=0; components<_meqn; components++)
                qOut[nodes*_meqn+components] = limititedValues[components];

        }
    }

    VecRestoreArrayRead(local_in, &u);
    VecRestoreArray(local_out, &v);

//    isInfinityOrNAN(local_in,"Limiter in NAN/INF");
//    isInfinityOrNAN(local_out,"Limiter in NAN/INF");

    DMLocalToGlobalBegin(_dm, local_out, INSERT_VALUES, q_limited);
    DMLocalToGlobalEnd(_dm, local_out, INSERT_VALUES, q_limited);

    DMRestoreLocalVector(_dm, &local_in);
    DMRestoreLocalVector(_dm, &local_out);

    VecDestroy(&local_out);
    VecDestroy(&local_in);
}


template<typename REAL>
void
WxTuAliabadiLimiter<REAL>::applyBc(int bcNum, REAL *xc, REAL *nx, REAL *q, REAL *qaux, REAL *AreaInts, REAL *qBC)
{
    // apply boundary conditions
    ApSubSolver<REAL>* ss = this->getParent()->getSubSolver( _bcSubSolvers.at(bcNum-1) );
    // cast this to the a grid BC and call step function
    dynamic_cast<WxGridBC<REAL>* >(ss)->applyToArray(xc,nx,q,qaux,AreaInts,qBC);
}

template <typename REAL>
PetscErrorCode
WxTuAliabadiLimiter<REAL>::isInfinityOrNAN(Vec f, std::string location)
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
//template class WxTuAliabadiLimiter<float>;
template class WxTuAliabadiLimiter<double>;
