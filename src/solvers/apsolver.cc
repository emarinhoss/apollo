// lib includes
#include <wxcreator.h>
#include <wxlogger.h>
#include <wxlogstream.h>
#include <wxtimer.h>
#include "petsc_compat.h"  // PETSc API compatibility for version 3.19+

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
//  PetscViewerDestroy(&_viewer);
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
//    PetscViewerCreate(PetscObjectComm((PetscObject)_dm), &_viewer);
//    PetscViewerSetType(_viewer, PETSCVIEWERVTK);
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

    unsigned nout;
    REAL tsize; // time between file output

    // write data to file before running main loop
    this->writeData(solution);

    // main solver loop
    tsize = _tend/_nout;

    for (unsigned i=0; i<_nout; ++i)
    {
        _tend_temp = (i+1)*tsize;
//        _dt_temp = tsize/(ceil(tsize/_dt));
        _dt_temp = _dt;
        tssolver->setTimeParameters(i*tsize, _tend_temp, _dt_temp);

        infStrm << "Advancing solution"
                << " from time " << (i*tsize)
                << " to " << _tend_temp
                << "..."
                << std::endl;

        WxTimer advTimer;
        // advance solution on each block by 'tsize'
        advTimer.startTimer();
        bool statusPetsc = tssolver->solve(solution);
        advTimer.stopTimer();

        // write solution to file
        this->writeData(solution);

        infStrm << "Advance completed in "
                << advTimer.timeElapsedAsString()
                << std::endl << std::endl;
    }
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
    tssolver = new WxPetscTimeSteppingSolver<REAL, ApSolver>(_dm, this, PETSC_COMM_WORLD);

    // DMCreateGlobalVector(_dm, &_usr.cg_vars);
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
//    PetscSF        sf;

    PetscFunctionBeginUser;
    MPI_Comm_rank(comm, &rank);
    PetscStrlen(_filename, &len);
    PetscStrncmp(&_filename[PetscMax(0,len-4)], extGmsh,   4, &isGmsh);
    PetscStrncmp(&_filename[PetscMax(0,len-4)], extExodus, 4, &isExodus);

    WxLogger *log = WxLogger::get("apollo-root.console");
    WxLogStream debStrm = log->getInfoStream();

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
//        PetscFinalize();
        exit(1);
    }

    // Distribute mesh over processes
    DM dmDist;
    DMPlexSetAdjacencyUseCone(*dm, PETSC_TRUE);
    DMPlexSetAdjacencyUseClosure(*dm, PETSC_FALSE);
    DMPlexDistribute(*dm, 1, NULL, &dmDist);
    if (dmDist){
        DMDestroy(dm);
        *dm   = dmDist;}
    // get any additional parameters for DM from the command line
    DMSetFromOptions(*dm);

    DM gdm;
    DMPlexConstructGhostCells(*dm, NULL, NULL, &gdm);
    DMDestroy(dm);
    *dm = gdm;

    // SetUp
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
    PetscFunctionBeginUser;
    PetscViewerCreate(PetscObjectComm((PetscObject)dm), viewer);
    PetscViewerSetType(*viewer, PETSCVIEWERVTK);
    PetscViewerFileSetName(*viewer, filename);
}

template<typename REAL>
void
ApSolver<REAL>::writeData(Vec X)
{
    PetscViewer viewer;
    std::stringstream ss; ss << _frameNum;
    std::string fname = this->runName() + "_" + ss.str() + ".vtu";
//    std::string fname = this->runName() + "_" + ss.str() + ".h5";
//    PetscViewerHDF5Open(PetscObjectComm((PetscObject)_dm),&fname[0],FILE_MODE_WRITE,&viewer);
    this->OutputVTK(_dm,&fname[0],&viewer);
    PetscBarrier((PetscObject) _dm);
    VecView(X,viewer);
    _frameNum += 1;
//    PetscBarrier((PetscObject) _dm);
    PetscViewerDestroy(&viewer);
}

template<typename REAL>
PetscErrorCode
ApSolver<REAL>::MonitorVTK(TS ts, PetscInt stepnum, PetscReal time, Vec X, void *ctx)
{
    WxLogger *log = WxLogger::get("apollo-root.console");
    WxLogStream infStrm = log->getInfoStream();

    // Adjust time-step
    if(fabs(_tend_temp-time)<_dt_temp){_dt_temp = fabs(_tend_temp-time);}
    REAL newdt = _dt_temp;
    MPI_Allreduce(&newdt, &_dt_temp, 1,
                  MPI_DOUBLE, MPI_MIN,
                  MPI_COMM_WORLD);
    TSSetTimeStep(ts,_dt_temp);
    PetscReal dtStep;
    TSGetTimeStep(ts,&dtStep);
//    PetscBarrier((PetscObject) _dm);
    infStrm << " Simulation dt = " << dtStep << " at t = " << time << std::endl;
    PetscFunctionReturn(0);
}

template<typename REAL>
PetscErrorCode
ApSolver<REAL>::ComputeRHSforTS(TS ts,PetscReal t,Vec global_in,Vec global_out,void *ctx)
{
    // Don't do anything if nothing to be done.
    if (_perStep.size() == 0) return 0;

    Vec X;
    VecDuplicate(global_in,&X);
    VecZeroEntries(global_out);
    WxStepperStatus<REAL> status;

    WxLogger *log = WxLogger::get("apollo-root.console");
    WxLogStream debStrm = log->getDebugStream();
    typename std::vector<ApSubSolverStep<REAL> >::iterator itr;

    for (itr=_perStep.begin(); itr!=_perStep.end(); ++itr)
    {
        PetscReal dtStep;
        TSGetTimeStep(ts,&dtStep);
        _dt_temp = dtStep;
        std::vector<std::string>::const_iterator ssitr;

        for (ssitr = itr->subSolvers.begin(); ssitr != itr->subSolvers.end(); ++ssitr)
        {
            debStrm << "  SubSolver " << *ssitr << std::endl;
            ApSubSolver<REAL> *ss = _subSolvers[*ssitr];
            // take this step
            status = ss->step(t,_dt_temp, global_in, X);
            //VecView(u,PETSC_VIEWER_STDOUT_WORLD);
            VecAXPY(global_out,1.0, X);
//            REAL newdt = status.getSuggestedDt();
            _dt_temp = fmin(status.getSuggestedDt(),_dt_temp);
        }
    }
    PetscBarrier((PetscObject) _dm);
    VecDestroy(&X);
    return 0;
}

// instantiations
template class ApSolver<float>;
template class ApSolver<double>;
