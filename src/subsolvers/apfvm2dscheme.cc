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
    VecDestroy(&qvars);
}

template <typename REAL>
void
ApFVM2Dscheme<REAL>::setup(const WxCryptSet& wxc)
{
  // call base class setup first
  ApSubSolver<REAL>::setup(wxc);

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
  // allocate memory for waves, speeds and fluctuations
  _s = alloc_1d<REAL>(_mwave); // wave speeds

  std::vector<WxAny> initArrays = wxc.template get<std::vector<WxAny> >("Initialize");
  std::vector<WxAny>::const_iterator iaitr;
  for (iaitr = initArrays.begin(); iaitr != initArrays.end(); ++iaitr)
    _initArrays.push_back(wx_any_cast<std::string>(*iaitr));

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

    std::vector<std::string>::const_iterator itr;
    for (itr = _initArrays.begin(); itr != _initArrays.end(); ++itr)
    {
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
            if(xc){xc[0] = d[0];}
        }
    }
    VecRestoreArray(out, &x);
}

template <typename REAL>
WxStepperStatus<REAL>
ApFVM2Dscheme<REAL>::step(REAL dt, Vec in, Vec out)
{

    WxStepperStatus<REAL> status;

    status.setStatus(true);
    status.setSuggestedDt(0.01);
    return status;

}

// instantiations
template class ApFVM2Dscheme<float>;
template class ApFVM2Dscheme<double>;
