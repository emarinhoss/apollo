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
ApDomainDecompCheck<REAL>::setup(const WxCryptSet& wxc, DM dm)
{
  ApSubSolver<REAL>::setup(wxc, dm);

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
ApDomainDecompCheck<REAL>::step(REAL dt, Vec in, Vec out)
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
