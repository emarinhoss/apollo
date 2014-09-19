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

  DMDestroy(&_dm);
  PetscViewerDestroy(&_viewer);
  TSDestroy(&_ts);
}

template <typename REAL>
void
ApSolver<REAL>::setStartFrame(unsigned frame)
{
  _startFrame = frame;
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
    _nout = wxc.template get<int>("Output_files");

    // time step to take
    _dt = wxc.template get<REAL>("Dt");

    // problem dimensions
    _dim = wxc.template get<int>("Dimensions");

    // grid to be used
    std::string fname = wxc.template get<std::string>("Gridname");
    _filename = &fname[0];

    // read and create mesh object
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
      ss->setup(sscs); // setup the subsolver
//      ss->declareTypes(); // get to declare its expected variable types
      _subSolvers.insert(SubSolverPair_t(sscs.name(), ss));
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
        sss.setup(ssscs);
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
        sss.setup(ssscs);
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
        sss.setup(ssscs);
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
        sss.setup(ssscs);
        _writeOnly.push_back(sss);
      }
    }

    // Set up data structure
    std::vector<WxAny> vars = wxc.template get<std::vector<WxAny> >("Variables");
    std::vector<WxAny> comp = wxc.template get<std::vector<WxAny> >("VariablesComponents");
    for(itr=vars.begin();itr!=vars.end();++itr)
    {
        std::string variable = wx_any_cast<std::string>(*itr);
        _fieldNames.push_back(variable);
    }
    for(itr=comp.begin();itr!=comp.end();++itr)
    {
        int num = wx_any_cast<int>(*itr);
        _fieldComponents.push_back(num);
    }
    if(_fieldComponents.size()!=_fieldNames.size())
    {
        std::cerr << "Variables vector and VariableComponents vector must be the same size." << std::endl;
        exit(1);
    }
    this->SetupLocalSpace(_dm,_usr);

    // Create viewer to output data into grid
    PetscViewerCreate(PetscObjectComm((PetscObject)_dm), &_viewer);
    PetscViewerSetType(_viewer, PETSCVIEWERVTK);

    // create time stepping scheme
    TSCreate(PETSC_COMM_WORLD, &_ts);
    TSSetType(_ts, TSSSP);
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
REAL
ApSolver<REAL>::getInitDt() const
{
  return _dt;
}

template <typename REAL>
void
ApSolver<REAL>::solve()
{

    // fetch stream for logging messages
    WxLogger *log = WxLogger::get("apollo-root.console");
    WxLogStream infStrm = log->getInfoStream();
    WxLogStream debStrm = log->getDebugStream();

    //REAL tcurr = this->getCurrentTime(); // starting time
    PetscScalar ftime;
    PetscInt nsteps;
    TSConvergedReason reason;

    // function that is to be used at every timestep
    // to display the iteration's progress.
    //TSMonitorSet(_ts,*MonitorVTK,NULL,NULL);

    TSSetDuration(_ts,_nout,_tend);
    TSSetInitialTimeStep(_ts,_tstart,_dt);
    TSSetFromOptions(_ts);
    TSSolve(_ts,_usr.solution);
    TSGetSolveTime(_ts,&ftime);
    TSGetTimeStepNumber(_ts,&nsteps);
    TSGetConvergedReason(_ts,&reason);

    infStrm << TSConvergedReasons[reason] << " at time "
            << ftime << " after "
            << nsteps << " steps"
            << std::endl;

}

template <typename REAL>
WxStepperStatus<REAL>
ApSolver<REAL>::step(REAL dt, Vec in, Vec out)
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
        WxStepperStatus<REAL> res = _subSolvers[*ssitr]->step(0.0, NULL, _usr.solution);
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

    // relate the timestepping scheme
    // with the data managenent object
    TSSetDM(_ts, _dm);
    TSSetRHSFunction(_ts,NULL,TSComputeRHSFunctionLinear,&_usr);

    // Solution vector
    DMCreateGlobalVector(_dm, &_usr.solution);
    PetscObjectSetName((PetscObject) _usr.solution, "solution");

    // set current time
    this->setCurrentTime(_tstart);
    // set frame number
    this->setStartFrame(0);

    std::string fname = this->runName() + "_0.vtu";
    this->setFilename_OutputVTK(fname);

    // initialize subsolvers
    typename SubSolverMap_t::iterator itr;
    for (itr = _subSolvers.begin(); itr != _subSolvers.end(); ++itr)
        itr->second->init();

    DMCreateGlobalVector(_dm, &_usr.cg_vars);
    // run startOnly subsolvers
    debStrm << "Running StartOnly steps...\n" << std::endl;
    startOnly();
}

template <typename REAL>
void
ApSolver<REAL>::advance(REAL tstart, REAL tend, REAL&dt)
{
    // don't do anything if nothing to do
    if (_perStep.size() == 0) return;

    WxLogger *log = WxLogger::get("apollo-root.console");
    WxLogStream debStrm = log->getDebugStream();
    WxLogStream infStrm = log->getInfoStream();

    REAL told, dtLast;
    REAL t = tstart, myDt = dt, dtNext, oldDt;
    unsigned nstep = 1, nsteps = 1;


    if (_useFixedDt)
      // the number of steps required to reach the output frame
      nsteps = unsigned(floor((tend-tstart)/dt));

    typename std::vector<ApSubSolverStep<REAL> >::iterator itr;

    //WxMsgBase& msg = this->getMsg();

    // loop advancing solution using adaptive time-stepping
    while (1)
    {
      dtNext = std::numeric_limits<REAL>::max();
      told = t;
      dtLast = myDt;

      // adjust dt to hit tend exactly if we are near the end of the
      // computation if fuzzy stepper is not used
      if (!_useFixedDt)
      {
        if (told+myDt>tend)
          myDt = tend-told;
      }

      redo:
      // ensure we have minimum time-step between processors
//      REAL myMinDt;
//      msg.allReduce(1, &myDt, &myMinDt, WX_MSG_MIN);
//      myDt = myMinDt;
      t = told + myDt;
      this->setCurrentTime(told);
      // advance solution by calling each substep
      for (itr=_perStep.begin(); itr!=_perStep.end(); ++itr)
      {
        REAL dtStep = itr->dtFrac*myDt;
        this->setDt(dtStep);

        std::vector<std::string>::const_iterator ssitr;
        for (ssitr = itr->subSolvers.begin(); ssitr != itr->subSolvers.end(); ++ssitr)
        {
          debStrm << " SubSolver " << *ssitr << std::endl;
          ApSubSolver<REAL> *ss = _subSolvers[*ssitr];

          // set time before calling step
          ss->setCurrentTime(told);
          ss->setDt(dtStep);
          // take this step
          WxStepperStatus<REAL> res = ss->step(dtStep, _usr.solution, _usr.solution);

          // check if step failed or succeeded.
          unsigned myStatus = res.getStatus();
          REAL newDt = res.getSuggestedDt();
          if (_useFixedDt)
            newDt = _dt;

          // unsigned status;
          // find if subsolver on processor failed
          // msg.allReduce(1, &myStatus, &status, WX_MSG_AND);
          if (myStatus)
          {
            dtNext = std::min(newDt, dtNext);
          }
          else
          {
            oldDt = myDt;
            myDt = newDt;

            infStrm << " **** Rejecting step " << nstep << ". Time step  " << oldDt << " too large"
                    << std::endl;
            if (_useFixedDt)
            {
              infStrm << "Please use " << myDt << ". Exiting simulation ... " << std::endl << std::endl;
              exit(1);
            }
            goto redo;
          }
        }
      }
      debStrm << " Step " << nstep  << " Time " << t  << " dt " << myDt << std::endl;

      // adjust time step to ensure we are getting proper time step
      myDt = dtNext;

      nstep += 1;
      // break if we are done
      if(_useFixedDt)
      {
        if(nstep == nsteps+1)
          // the fuzzy stepper break is determined based on how many
          // steps have been taken between tstart and tend, and not
          // on how close the current time with tend. This is so
          // because it has been assumed that the time step will
          // always be constant throughout the simulation
          break;
      }
      else
      {
        if( (tend-t) < 5*tend*std::numeric_limits<REAL>::epsilon())
          break;
      }
    }
    dt = dtLast; // copy last time step
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
    DMPlexDistribute(*dm, "metis", 1, NULL, &dmDist);
    if (dmDist){
        DMDestroy(dm);
        *dm   = dmDist;}
    // get any additional parameters for DM from the command line
    DMSetFromOptions(*dm);

    DM gdm;
    DMPlexConstructGhostCells(*dm, NULL, NULL, &gdm);
    DMDestroy(dm);
    *dm = gdm;

    // SetUp the data structures
    DMSetUp(*dm);
    PetscObjectSetName((PetscObject) *dm, "Mesh");
}

template<typename REAL>
void
ApSolver<REAL>::SetupLocalSpace(DM dm, UserContext usr)
{
    PetscSection   stateSection;
    PetscInt       cStart, cEnd, c;

    DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd);
    PetscSectionCreate(PetscObjectComm((PetscObject)dm), &stateSection);    
    PetscSectionSetNumFields(stateSection,_fieldNames.size());
    for(int k=0; k<_fieldNames.size(); ++k)
    {
        PetscSectionSetFieldComponents(stateSection,k,_fieldComponents[k]);
        std::string tt=_fieldNames[k]; char *fieldN = &tt[0];
        PetscSectionSetFieldName(stateSection,0,fieldN);
    }

    PetscSectionSetChart(stateSection, cStart, cEnd);

    for (c = cStart; c < cEnd; ++c)
    {
        for(int i=0; i<_fieldNames.size(); i++)
        {
            PetscSectionSetFieldDof(stateSection,c,i,_fieldComponents[i]);
            PetscSectionSetDof(stateSection, c, _fieldComponents[i]);
        }
    }

    PetscSectionSetUp(stateSection);
    DMSetDefaultSection(dm,stateSection);
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
PetscErrorCode ApSolver<REAL>::MonitorVTK(TS ts, PetscInt stepnum, PetscReal time, Vec X, void *ctx)
{
    PetscViewer viewer;

    if ((stepnum == -1) ^ (stepnum % _nout == 0))
    {
        if(stepnum == -1) {/* Final time is not multiple of normal time interval, write it anyway */
          TSGetTimeStepNumber(ts,&stepnum);}

        std::stringstream ss; ss << stepnum;
        std::string fname = this->runName() + "_" + ss.str() + ".vtu";
        OutputVTK(_dm,&fname[0],&viewer);
        VecView(X,viewer);
        PetscViewerDestroy(&viewer);
      }

    PetscFunctionReturn(0);
}

// instantiations
template class ApSolver<float>;
template class ApSolver<double>;
