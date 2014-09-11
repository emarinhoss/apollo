// Apollo includes
#include "apdomaindecompcheck.h"

// WarpX lib includes
#include <wxcreator.h>
#include <wxmath.h>

// std includes
#include <string>
#include <iostream>

template <typename REAL>
void
ApDomainDecompCheck<REAL>::setup(const WxCryptSet& wxc)
{
  ApSubSolver<REAL>::setup(wxc);

//  _q = wxc.template get<REAL>("charge");
//  _mi = wxc.template get<REAL>("ionmass");
//  _me = wxc.template get<REAL>("elcmass");
//  _eps0 = wxc.template get<REAL>("epsilon0");

  // check if this simulation has transport
//  if (wxc.has("isTransport"))
//  {
//    std::string flg = wxc.template get<std::string>("isTransport");
//    if (flg == "true")
//    {
//      _transport=true;
//      _gas_gamma = wxc.template get<REAL>("gas_gamma");
//      _k = wxc.template get<REAL>("boltz");
//      _dx = wxc.template get<REAL>("dx");
//    }
//  }

}

template <typename REAL>
void
ApDomainDecompCheck<REAL>::init()
{
    _dm  = this->getParent()->getdatamanagment();
    _usr = this->getParent()->getusercontext();

}

template <typename REAL>
WxStepperStatus<REAL>
ApDomainDecompCheck<REAL>::step(REAL dt)
{
    PetscInt cStart, cEnd, rank;
    PetscScalar *x;
    Vec X;
    MPI_Comm comm;
    PetscViewer viewer = this->getParent()->getViewer();

    DMGetGlobalVector(_dm, &X);
    PetscObjectSetName((PetscObject) X, "partition");
    VecGetArray(X, &x);
    DMPlexGetHeightStratum(_dm, 0, &cStart, &cEnd);
    PetscObjectGetComm((PetscObject)_dm,&comm);
    MPI_Comm_rank(comm, &rank);

    for(int c = cStart; c < cEnd; ++c)
    {
        PetscScalar *xc;
        DMPlexPointGlobalRef(_dm,c,x,&xc);
        if (xc){xc[0] = rank;}
    }

    VecRestoreArray(X, &x);
    VecView(X,viewer);
    return true;

}

// instantiations
template class ApDomainDecompCheck<float>;
template class ApDomainDecompCheck<double>;
