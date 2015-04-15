// lib includes
#include <wxcreator.h>
#include <wxlogger.h>
#include <wxlogstream.h>
#include <wxtimer.h>

// solver includes
#include "apsolver.h"

// std includes
#include <cmath>
#include <ctime>
#include <iostream>
#include <string>
#include <vector>
#include <limits>


template <typename REAL>
ApSolver<REAL>::ApSolver(const std::string& name)
  : WxSolverBase<REAL>(name)
{
}

template <typename REAL>
ApSolver<REAL>::~ApSolver()
{
  // delete subsolvers
  typename SubSolverMap_t::iterator ssItr;
  for (ssItr=_subSolvers.begin(); ssItr!=_subSolvers.end(); ++ssItr)
    delete ssItr->second;

  delete tssolver;
//  DMDestroy(&_dm);
  VecDestroy(&solution);
  PetscViewerDestroy(&_viewer);
}

template <typename REAL>
void
ApSolver<REAL>::setup(const WxCryptSet& wxc)
{
    WxLogStream debStrm = WxLogger::get("apollo-root.console")->getDebugStream();

    // setup our parent first
    WxSolverBase<REAL>::setup(wxc);
    // start and end time for simulation
    std::vector<WxAny> times = wxc.template get<std::vector<WxAny> >("Time");
    _tstart = wx_any_cast<REAL>(times[0]);
    _tend = wx_any_cast<REAL>(times[1]);

    // number of files to write
     _nout = wxc.template get<int>("Output_files"); _usr.nout = _nout;

    // time step to take
    _dt = wxc.template get<REAL>("Dt");

    // problem dimensions
    _dim = wxc.template get<int>("Dimensions");

    // total number of fields
    _fieldsNum = wxc.template get<int>("NumFields");

    // grid to be used
    std::string fname = wxc.template get<std::string>("Gridname");
    _filename = &fname[0];

    // read and create mesh objectwxpnodaldgfunctions
    this->createMesh(PETSC_COMM_WORLD,&_dm);

    std::vector<std::string>::const_iterator i;

    // initialize all subsolvers
    std::vector<std::string> subsolverCS = wxc.getNamesOfType("ApSubSolver");
    for (i=subsolverCS.begin(); i!=subsolverCS.end(); ++i)
    {
      const WxCryptSet& sscs = wxc.getSet(*i);
      std::string kind = sscs.get<std::string>("Kind");
      ApSubSolver<REAL> *ss = WxCreatorMap<ApSubSolver<REAL> >::getNew(kind);
      ss->setIo(this->getIo());
      ss->setMsg(this->getMsg());
      ss->setParent(this);
      ss->setup(sscs,_dm); // setup the subsolver
      _subSolvers.insert(SubSolverPair_t(sscs.name(), ss));
      if (sscs.has("DataStructure"))
      {
          std::vector<WxAny> datastruc = sscs.get<std::vector<WxAny> >("DataStructure");
          _fieldsNumber.push_back(wx_any_cast<int>(datastruc[0]));
          _fieldsName.push_back(wx_any_cast<std::string>(datastruc[1]));
          _fieldsComponents.push_back(wx_any_cast<int>(datastruc[2]));
      }
    }

    // initialize various solver sequences
    const WxCryptSet sseqcs = wxc.getSet("SolverSequence");
    std::vector<WxAny>::const_iterator itr;
    // list of StartOnly subsolver steps
    if (sseqcs.has("StartOnly"))
    {
      std::vector<WxAny> list = sseqcs.template get<std::vector<WxAny> >("StartOnly");
      for (itr=list.begin(); itr!=list.end(); ++itr)
      {
        std::string name = wx_any_cast<std::string>(*itr);
        const WxCryptSet& ssscs = wxc.getSet(name);
        ApSubSolverStep<REAL> sss;
        sss.setup(ssscs,_dm);
        _startOnly.push_back(sss);
      }
    }
    // list of EndOnly subsolver steps
    if (sseqcs.has("EndOnly"))
    {
      std::vector<WxAny> list = sseqcs.template get<std::vector<WxAny> >("EndOnly");
      for (itr=list.begin(); itr!=list.end(); ++itr)
      {
        std::string name = wx_any_cast<std::string>(*itr);
        const WxCryptSet& ssscs = wxc.getSet(name);
        ApSubSolverStep<REAL> sss;
        sss.setup(ssscs,_dm);
        _endOnly.push_back(sss);
      }
    }
    // sequence of steps at each time step
    if (sseqcs.has("PerStep"))
    {
      std::vector<WxAny> list = sseqcs.template get<std::vector<WxAny> >("PerStep");
      for (itr=list.begin(); itr!=list.end(); ++itr)
      {
        std::string name = wx_any_cast<std::string>(*itr);
        const WxCryptSet& ssscs = wxc.getSet(name);
        ApSubSolverStep<REAL> sss;
        sss.setup(ssscs,_dm);
        _perStep.push_back(sss);
      }
    }
    // sequence of steps at before writing out data
    if (sseqcs.has("WriteOnly"))
    {
      std::vector<WxAny> list = sseqcs.template get<std::vector<WxAny> >("WriteOnly");
      for (itr=list.begin(); itr!=list.end(); ++itr)
      {
        std::string name = wx_any_cast<std::string>(*itr);
        const WxCryptSet& ssscs = wxc.getSet(name);
        ApSubSolverStep<REAL> sss;
        sss.setup(ssscs,_dm);
        _writeOnly.push_back(sss);
      }
    }

    // Setup data structure
    this->SetupLocalSpace(&_dm);

    // Create viewer to output data into grid
    PetscViewerCreate(PetscObjectComm((PetscObject)_dm), &_viewer);
    PetscViewerSetType(_viewer, PETSCVIEWERVTK);
}

template <typename REAL>
ApSubSolver<REAL>*
ApSolver<REAL>::getSubSolver(const std::string& name)
{
  typename SubSolverMap_t::iterator i;
  i = _subSolvers.find(name);
  if (i != _subSolvers.end())
    return i->second;
  WxExcept wxe;
  wxe << "WxSolver::getSubSolver : SubSolver " << name << " not found";
  throw wxe;
}


template <typename REAL>
void
ApSolver<REAL>::solve()
{

    // fetch stream for logging messages
    WxLogger *log = WxLogger::get("apollo-root.console");
    WxLogStream infStrm = log->getInfoStream();

    // solve equation system
    bool statusPetsc = tssolver->solve(solution);

}

template <typename REAL>
WxStepperStatus<REAL>
ApSolver<REAL>::step(REAL t, REAL dt, Vec in, Vec out)
{
  return WxStepperStatus<REAL>();
}

template <typename REAL>
void
ApSolver<REAL>::startOnly()
{
    WxLogger *log = WxLogger::get("apollo-root.console");
    WxLogStream debStrm = log->getDebugStream();

    typename std::vector<ApSubSolverStep<REAL> >::iterator itr;
    for (itr = _startOnly.begin(); itr!=_startOnly.end(); ++itr)
    {
      // run the subsolver step
      std::vector<std::string>::const_iterator ssitr;
      for (ssitr = itr->subSolvers.begin(); ssitr != itr->subSolvers.end(); ++ssitr)
      {
        debStrm << " SubSolver " << *ssitr << std::endl;

        _subSolvers[*ssitr]->setCurrentTime(this->getCurrentTime());
        _subSolvers[*ssitr]->setDt(0.0);
        WxStepperStatus<REAL> res = _subSolvers[*ssitr]->step(0.0, 0.0, NULL, solution);
        if (res.getStatus() == false)
        {
          WxExcept wxe("Subsolver ");
          wxe << *ssitr << " failed" << std::endl;
          throw wxe;
        }
      }
    }
}

template <typename REAL>
void
ApSolver<REAL>::init()
{
    WxLogger *log = WxLogger::get("apollo-root.console");
    WxLogStream debStrm = log->getDebugStream();
    REAL suggestedDt;

    // Solution vector
    DMCreateGlobalVector(_dm, &solution);
    PetscObjectSetName((PetscObject) solution, "solution");

    // time step
    // REAL dtInit = 1.e6;
    _frameNum = 0;

    // initialize subsolvers
    typename SubSolverMap_t::iterator itr;
    for (itr = _subSolvers.begin(); itr != _subSolvers.end(); ++itr)
        itr->second->init(_dt, solution);

    typename std::vector<ApSubSolverStep<REAL> >::iterator itrr;
    for (itrr = _perStep.begin(); itrr!=_perStep.end(); ++itrr)
    {
      // run the subsolver step
      std::vector<std::string>::const_iterator ssitr;
      for (ssitr = itrr->subSolvers.begin(); ssitr != itrr->subSolvers.end(); ++ssitr)
      {
        suggestedDt = _subSolvers[*ssitr]->getDt();
        _dt = fmin(_dt,suggestedDt);
      }
    }

    // Initialize the timestepping solver
    tssolver = new WxPetscTimeSteppingSolver<REAL, ApSolver>(_dm, this, PetscObjectComm((PetscObject)_dm),_tstart,_tend,_dt);

    DMCreateGlobalVector(_dm, &_usr.cg_vars);
    // run startOnly subsolvers
    debStrm << "Running StartOnly steps...\n" << std::endl;
    startOnly();
}

template <typename REAL>
void
ApSolver<REAL>::createMesh(MPI_Comm comm, DM *dm)
{
    const char    *extGmsh     = ".msh";
    const char    *extExodus   = ".exo";
    size_t         len;
    PetscBool      isGmsh, isExodus;
    PetscMPIInt    rank;

    PetscFunctionBeginUser;
    MPI_Comm_rank(comm, &rank);
    PetscStrlen(_filename, &len);
    PetscStrncmp(&_filename[PetscMax(0,len-4)], extGmsh,   4, &isGmsh);
    PetscStrncmp(&_filename[PetscMax(0,len-4)], extExodus, 4, &isExodus);

    WxLogger *log = WxLogger::get("apollo-root.console");
    WxLogStream debStrm = log->getDebugStream();

    debStrm << "Reading grid file --> " << _filename << std::endl;

    if (isGmsh)
    {
        PetscViewer viewer;

        PetscViewerCreate(comm, &viewer);
        PetscViewerSetType(viewer, PETSCVIEWERASCII);
        PetscViewerFileSetMode(viewer, FILE_MODE_READ);
        PetscViewerFileSetName(viewer, _filename);
        DMPlexCreateGmsh(comm, viewer, PETSC_TRUE, dm);
        PetscViewerDestroy(&viewer);
    }
    else if (isExodus)
    {
        DMPlexCreateExodusFromFile(comm, _filename, PETSC_TRUE, dm);
    }
    else
    {
        std::cerr << "Mesh input filename " << _filename << " not found.  Exiting." << std::endl;
        exit(1);
    }

    // Distribute mesh over processes
    DM dmDist;
    debStrm << "Partitioning the domain using --> Metis "  << std::endl;
    //debStrm << "Partitioning the domain using --> " << _partitioner  << std::endl;
    DMPlexDistribute(*dm,"metis", 0, NULL, &dmDist);
    if (dmDist){
        DMDestroy(dm);
        *dm   = dmDist;}
    // get any additional parameters for DM from the command line
    DMSetFromOptions(*dm);

//    DM gdm;
//    DMPlexConstructGhostCells(*dm, NULL, NULL, &gdm);
//    DMDestroy(dm);
//    *dm = gdm;

    // SetUp the data structures
    DMSetUp(*dm);
    PetscObjectSetName((PetscObject) *dm, "Mesh");
}

template<typename REAL>
void
ApSolver<REAL>::SetupLocalSpace(DM *dm)
{
    PetscSection   stateSection;
    PetscInt       cStart, cEnd, c;

    DMPlexGetHeightStratum(*dm, 0, &cStart, &cEnd);
    PetscSectionCreate(PetscObjectComm((PetscObject)*dm), &stateSection);
    PetscSectionSetNumFields(stateSection,_fieldsNum);

    for(unsigned k=0; k<_fieldsNum; k++)
    {
        PetscSectionSetFieldComponents(stateSection,_fieldsNumber[k],_fieldsComponents[k]);
        std::string name = _fieldsName[k];
        PetscSectionSetFieldName(stateSection,_fieldsNumber[k],&name[0]);
    }
    PetscSectionSetChart(stateSection, cStart, cEnd);

    for (c = cStart; c < cEnd; ++c)
    {
        for(unsigned kk=0; kk<_fieldsNum; kk++)
        {
            PetscSectionSetFieldDof(stateSection,c,_fieldsNumber[kk],_fieldsComponents[kk]);
            PetscSectionSetDof(stateSection, c, _fieldsComponents[kk]);
        }
    }

    PetscSectionSetUp(stateSection);
    DMSetDefaultSection(*dm,stateSection);
    PetscSectionDestroy(&stateSection);
}

template<typename REAL>
void
ApSolver<REAL>::OutputVTK(DM dm, char *filename, PetscViewer *viewer)
{
    PetscViewerCreate(PetscObjectComm((PetscObject)dm), viewer);
    PetscViewerSetType(*viewer, PETSCVIEWERVTK);
    PetscViewerFileSetName(*viewer, filename);
}

template<typename REAL>
PetscErrorCode
ApSolver<REAL>::MonitorVTK(TS ts, PetscInt stepnum, PetscReal time, Vec X, void *ctx)
{
    WxLogger *log = WxLogger::get("apollo-root.console");
    WxLogStream infStrm = log->getInfoStream();

    if ((stepnum == -1) ^ (stepnum % _nout == 0))
    {
        PetscViewer viewer;
        if(stepnum == -1) {/* Final time is not multiple of normal time interval, write it anyway */
          TSGetTimeStepNumber(ts,&stepnum);}

        std::stringstream ss; ss << _frameNum;
        std::string fname = this->runName() + "_" + ss.str() + ".vtu";
        //PetscViewerHDF5Open(PetscObjectComm((PetscObject)ts),&fname[0],FILE_MODE_WRITE,&viewer);
        this->OutputVTK(_dm,&fname[0],&viewer);
        VecView(X,viewer);
        _frameNum += 1;
        PetscViewerDestroy(&viewer);
      }

    // Adjust time-step
    if(fabs(_tend-time)<_dt){_dt = fabs(_tend-time);}
    TSSetTimeStep(ts,_dt);
    PetscReal dtStep;
    TSGetTimeStep(ts,&dtStep);
    infStrm << " Current simulation time-step is " << dtStep << " and current time is " << time << std::endl;
    PetscFunctionReturn(0);
}

template<typename REAL>
PetscErrorCode
ApSolver<REAL>::ComputeRHSforTS(TS ts,PetscReal t,Vec u,Vec F,void *ctx)
{
    // don't do anything if nothing to do
    if (_perStep.size() == 0) return 0;

    Vec X;
    VecDuplicate(u,&X);
    VecZeroEntries(F);
    WxStepperStatus<REAL> status;

    WxLogger *log = WxLogger::get("apollo-root.console");
    WxLogStream debStrm = log->getDebugStream();
    WxLogStream infStrm = log->getInfoStream();
    typename std::vector<ApSubSolverStep<REAL> >::iterator itr;

    for (itr=_perStep.begin(); itr!=_perStep.end(); ++itr)
    {
        PetscReal dtStep;
        TSGetTimeStep(ts,&dtStep);
        _dt = dtStep;
        std::vector<std::string>::const_iterator ssitr;

        for (ssitr = itr->subSolvers.begin(); ssitr != itr->subSolvers.end(); ++ssitr)
        {
            debStrm << "  SubSolver " << *ssitr << std::endl;
            ApSubSolver<REAL> *ss = _subSolvers[*ssitr];
            // take this step
            status = ss->step(t,_dt, u, X);
            //VecView(u,PETSC_VIEWER_STDOUT_WORLD);
            VecAXPY(F,1.0, X);
            _dt = fmin(status.getSuggestedDt(),_dt);
        }
    }
    VecDestroy(&X);
    return 0;
}

// instantiations
template class ApSolver<float>;
template class ApSolver<double>;
