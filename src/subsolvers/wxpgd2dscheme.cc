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
    delete [] _ql;
    delete [] _qr;
    delete [] _fl;
    delete [] _fr;
    delete [] _src;
    delete [] _df;
}

template <typename REAL>
void
WxpDG2Dscheme<REAL>::setup(const WxCryptSet& wxc, DM dm)
{
  // call base class setup first
  ApSubSolver<REAL>::setup(wxc, dm);

  _spatialOrder = wxc.template get<int>("spatialOrder");

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

  // allocate memory for waves, speeds and fluctuations
  _s = alloc_1d<REAL>(_mwave); // wave speeds

  _dataStruct = wxc.template get<std::vector<WxAny> >("DataStructure");
  // add number of components
  _dataStruct.push_back(_meqn);
  // add total number of dofs
  _dataStruct.push_back((_spatialOrder+1)*(_spatialOrder+2)/2);

  // create function pointer for initial condition
  const WxCryptSet& initCS = wxc.getSet("InitialCondition");
  std::string kind;
  kind = initCS.template get<std::string>("Kind");
  _initFunc = WxCreatorMap<WxFunction<REAL> >::getNew(kind);
  // setup this function
  _initFunc->setup(initCS);

  _quad = new WxpDGGeometry<REAL>(dm, _meqn, _spatialOrder);
}

template <typename REAL>
void
WxpDG2Dscheme<REAL>::init(Vec out)
{
    DM dm;
    VecGetDM(out, &dm);
    PetscSection stateSection;
    DMGetDefaultSection(dm, &stateSection);
    PetscInt secComp = wx_any_cast<int>(_dataStruct[2]);

    // Finite Volume geometry variables
    PetscReal centroid[3], normal[3], vol;
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
WxpDG2Dscheme<REAL>::step(REAL dt, Vec in, Vec out)
{

}

// instantiations
template class WxpDG2Dscheme<float>;
template class WxpDG2Dscheme<double>;
