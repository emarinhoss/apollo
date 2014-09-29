#include "apfvm2dscheme.h"

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
ApFVM2Dscheme<REAL>::~ApFVM2Dscheme() {
    delete [] _ql;
    delete [] _qr;
    delete [] _fl;
    delete [] _fr;
    delete [] _src;
    delete [] _df;
    delete [] _s;
    delete [] _qauxl;
    delete [] _qauxr;
    delete [] _amdqx;
    delete [] _apdqx;
    delete [] _amdq;
    delete [] _apdq;
    free_2d_c(_wave, _meqn, _meqn);
    free_2d_c(_waveax, _meqn, _meqn);
}

template <typename REAL>
void
ApFVM2Dscheme<REAL>::setup(const WxCryptSet& wxc, DM dm)
{
  // call base class setup first
  ApSubSolver<REAL>::setup(wxc, dm);

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

  _ql = alloc_1d<REAL>(_meqn);
  _qr = alloc_1d<REAL>(_meqn);
  _qauxl = alloc_1d<REAL>(_meqn);
  _qauxr = alloc_1d<REAL>(_meqn);
  _df = alloc_1d<REAL>(_meqn); // jump
  _src = alloc_1d<REAL>(_meqn);
  _fl = alloc_1d<REAL>(_meqn);
  _fr = alloc_1d<REAL>(_meqn);
  _wave = alloc_2d_c<REAL>(_meqn, _mwave); // waves
  _waveax = alloc_2d_c<REAL>(_meqn, _mwave); // waves
  _amdqx = alloc_1d<REAL>(_meqn);
  _apdqx = alloc_1d<REAL>(_meqn);
  _amdq = alloc_1d<REAL>(_meqn);
  _apdq = alloc_1d<REAL>(_meqn);
  // allocate memory for waves, speeds and fluctuations
  _s = alloc_1d<REAL>(_mwave); // wave speeds

  _dataStruct = wxc.template get<std::vector<WxAny> >("DataStructure");

  // create function pointer for initial condition
  const WxCryptSet& initCS = wxc.getSet("InitialCondition");
  std::string kind;
  kind = initCS.template get<std::string>("Kind");
  _initFunc = WxCreatorMap<WxFunction<REAL> >::getNew(kind);
  // setup this function
  _initFunc->setup(initCS);


}

template <typename REAL>
void
ApFVM2Dscheme<REAL>::init(Vec out)
{
    DM dm;
    VecGetDM(out, &dm);
    PetscSection stateSection;
    DMGetDefaultSection(dm, &stateSection);
    PetscInt secComp = wx_any_cast<int>(_dataStruct[2]);

    // Finite Volume geometry variables
    PetscReal centroid[3], normal[3], vol;
    PetscScalar *x;

    PetscInt cStart, cEnd, c, cEndInterior;

    // Independent variables (t,x,y,z)
    REAL txo[5];
    txo[0] = this->getCurrentTime();

    // results returned by the initialization function
    REAL *d = new REAL[1]; // TODO: change 1 to the number of components

    // ****
    // Initialize the solution vector
    // ****
    // Get cells in this processor
    DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd);
    DMPlexGetHybridBounds(dm, &cEndInterior, NULL, NULL, NULL);
    VecGetArray(out, &x);

    for (c = cStart; c < cEndInterior; ++c)
    {
        DMPlexComputeCellGeometryFVM(dm, c, &vol, centroid, normal);
        txo[1] = centroid[0];
        txo[2] = centroid[1];
        txo[3] = centroid[2];

        _initFunc->func(3, txo, d);
        PetscScalar *xc;

        // reference this cell to the proper location on the solution
        // vector
        DMPlexPointGlobalRef(dm,c,x,&xc);
        // assign value returned by the initialization function
        // to the solution vector
        if(xc){
            for(unsigned kk=0; kk<secComp; kk++)
                xc[kk] = d[kk];
        }
    }
    VecRestoreArray(out, &x);

}

template <typename REAL>
WxStepperStatus<REAL>
ApFVM2Dscheme<REAL>::step(REAL dt, Vec in, Vec out)
{
//    VecView(in,PETSC_VIEWER_STDOUT_WORLD);
    Vec locX, locCorr;
    WxStepperStatus<REAL> status;
    const PetscScalar *x;
    PetscScalar *ot;
    // data for right/left rotated data
    std::vector<REAL> qrLocal(_meqn), qlLocal(_meqn);
    // data for wave in local and global coordinates
    std::vector<REAL> waveLocal(_meqn), waveGlobal(_meqn);

    // get the data managent object from the input vector
    VecGetDM(in,&_dm);

    // ghost cell label
    DMLabel ghostLabel, faceS;
    DMPlexGetLabel(_dm, "ghost", &ghostLabel);
    DMPlexGetLabel(_dm, "Face Sets", &faceS);

    // create local vector
    DMGetLocalVector(_dm, &locX);
    //DMGetLocalVector(_dm, &locCorr);

    // zero entries of the vectors that will be used to store
    // information
    VecZeroEntries(locX);
    //VecZeroEntries(locCorr);
    VecZeroEntries(out);

    // get local values of the global vector in into locX
    DMGlobalToLocalBegin(_dm, in, INSERT_VALUES, locX);
    DMGlobalToLocalEnd(_dm, in, INSERT_VALUES, locX);
//    VecView(locX,PETSC_VIEWER_STDOUT_WORLD);

    // get the field number from the data strucure
    int fieldnum = wx_any_cast<int>(_dataStruct[0]);
    DMGetField(_dm, fieldnum, (PetscObject *) &_fvm);
    PetscFVGetLimiter(_fvm, &_lim);

    // get start and end of faces
    PetscInt fStart, fEnd;
    DMPlexGetHeightStratum(_dm, 1, &fStart, &fEnd);
    //VecGetDM(_fvgeom->facegeom, &dmFace);
    //VecGetDM(_fvgeom->cellgeom, &dmCell);
    //VecGetArrayRead(_fvgeom->facegeom, &facegeom);
    //VecGetArrayRead(_fvgeom->cellgeom, &cellgeom);
    VecGetArrayRead(locX, &x);
    //VecGetArray(locCorr, &f);
    VecGetArray(out, &ot);

    for(PetscInt face = fStart; face < fEnd; ++face)
    {
        // is this a ghost cell face?
        PetscInt ghost;
        DMLabelGetValue(faceS, face, &ghost);

        if(ghost>=0)
        {
            // NOT QUITE SURE WHY THIS IS NEED ---- DOUBLE CHECK
            PetscBool boundary;
            DMPlexIsBoundaryPoint(_dm, face, &boundary);
            if (!boundary)
            {
                // Get the cells that support this face
                const PetscInt *cells;
                //CellGeom *cgL, *cgR;
                PetscReal cgL[3], cgR[3];
                DMPlexGetSupport(_dm, face, &cells);
                PetscReal volumeL, volumeR;
                DMPlexComputeCellGeometryFVM(_dm, cells[0], &volumeL, cgL, NULL);
                DMPlexComputeCellGeometryFVM(_dm, cells[1], &volumeR, cgR, NULL);

                // get the geometry of this face
                //FaceGeom *fg;
                PetscReal faceNormal[3];
                PetscReal area;
                DMPlexComputeCellGeometryFVM(_dm, face, &area, NULL, faceNormal);
                REAL normal[3]; for(unsigned k=0; k<3; k++){normal[k]=faceNormal[k];}

                // read conserved variable to the right and left
                // of this face
                REAL *qL, *qR;
                DMPlexPointLocalRead(_dm, cells[0], x, &qL);
                DMPlexPointLocalRead(_dm, cells[1], x, &qR);

                // get the centroid coordinates
                REAL xl[4], xr[4];
                for(unsigned k=0; k<3; k++){
                    xl[k] = cgL[k];
                    xr[k] = cgR[k];
                }

                _eqnSet.rotateToLocalFrame(normal, qL, &qrLocal[0]);
                _eqnSet.rotateToLocalFrame(normal, qR, &qlLocal[0]);

                // compute jump
                for (unsigned m=0; m<_meqn; ++m)
                  _df[m] = qrLocal[m] - qlLocal[m]; // jump in q across interface

                // call Riemann problem solver to get waves
                _eqnSet.riemann(0, xl, xr, &qlLocal[0], &qrLocal[0],
                                normal, 0, _df, _wave, _s, _amdqx, _apdqx);

                // rotate wave back to global coordinate
                for (unsigned mw=0; mw<_mwave; ++mw)
                {
                  // copy wave into temporary array
                  for (unsigned m=0; m<_meqn; ++m)
                    waveLocal[m] = _wave[m][mw];
                  // rotate it
                  _eqnSet.rotateToGlobalFrame(normal, &waveLocal[0], &waveGlobal[0]);
                  // copy it into array
                  for (unsigned m=0; m<_meqn; ++m)
                    _waveax[m][mw] = waveGlobal[m];
                }

                // now compute fluctuations in global coordinates
                for (unsigned m=0; m<_meqn; ++m)
                {
                  _amdq[m] = 0.0;
                  _apdq[m] = 0.0;
                  for (unsigned mw=0; mw<_mwave; ++mw)
                  {
                    if (_s[mw] < 0.0)
                      _amdq[m] += _s[mw]*_waveax[m][mw];
                    else
                      _apdq[m] += _s[mw]*_waveax[m][mw];
                  }
                }

                // compute Area/Volume ratios ratios
                REAL areaVol = area/volumeL;
                REAL areaVol1 = area/volumeR;

                // write first order Gudonov updates
                PetscScalar *fL, *fR;
                DMPlexPointGlobalRef(_dm, cells[0], ot, &fL);
                DMPlexPointGlobalRef(_dm, cells[1], ot, &fR);
                for(unsigned m=0; m<_meqn; ++m)
                {
                    fL[m] += areaVol*_apdq[m];
                    fR[m] += areaVol1*_amdq[m];
                }

//                /** LIMITER GOES SOMEWHERE HERE */

//                // compute second order corrections to fluxes
//                REAL cellVol = 0.5*(volumeL+volumeR); // average volume
//                REAL areaVavg = area/cellVol;


//                for (unsigned m=0; m<_meqn; ++m)
//                {
//                  _fsx[m] = 0.0;
//                  for (unsigned mw=0; mw<_mwave; ++mw)
//                  {
//                    REAL sabs = fabs(_s[mw]);
//                    REAL corr = 0.5*sabs*(1.0 - sabs*dt*areaVavg)*_waveax[m][mw];
//                    _fsx[m] += corr; // compute second order correction
//                  }
//                }
            }
        }
    }

    PetscInt eStart, eEnd, eEndInterior;
    DMPlexGetHeightStratum(_dm, 0, &eStart, &eEnd);
    DMPlexGetHybridBounds(_dm, &eEndInterior, NULL, NULL, NULL);
    for(unsigned kk=eEndInterior; kk<eEnd; kk++)
    {
        PetscScalar *oo;
        DMPlexPointGlobalRef(_dm, kk, ot, &oo);
        oo[0] = 0.0;
    }


    DMRestoreLocalVector(_dm, &locX);
    VecRestoreArray(out, &ot);
    //VecView(out,PETSC_VIEWER_STDOUT_WORLD);

    status.setStatus(true);
    status.setSuggestedDt(0.01);
    return status;

}

// instantiations
template class ApFVM2Dscheme<float>;
template class ApFVM2Dscheme<double>;
