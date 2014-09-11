// lib includes
#include <wxcreator.h>
#include <wxlogger.h>
#include <wxlogstream.h>
#include <wxtimer.h>

// solver includes
#include "apsolver.h"
//#include <wxwriteonlysubsolver.h>

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

    // grid partinioner to be used
    _partitioner = "metis";
    if (wxc.has("Partinioner_name"))
    {
        std::string pname = wxc.template get<std::string>("Partinioner_name");
        _partitioner = &pname[0];
    }

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

    this->SetupLocalSpace(_dm,_usr);
    PetscViewerCreate(PetscObjectComm((PetscObject)_dm), &_viewer);
    PetscViewerSetType(_viewer, PETSCVIEWERVTK);
    //PetscViewerFileSetName(_viewer, "test.vtu");
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

  REAL tcurr = this->getCurrentTime(); // starting time
  unsigned frame = _startFrame; // frame number

  //REAL temp_var;

  //unsigned nout;
  //REAL tsize; // time between file output
  // write data to file before running main loop
  //this->OutputVTK(_dm, frame);
  //PetscViewerFileSetName(_viewer, "test.vtu");


  // main solver loop

  frame += 1;

}

template <typename REAL>
WxStepperStatus<REAL>
ApSolver<REAL>::step(REAL dt)
{
  return WxStepperStatus<REAL>();
}

template <typename REAL>
void
ApSolver<REAL>::OutputVTK(DM dm, unsigned frame)
{
    // create new file
    PetscViewerFileSetName(_viewer, "test.vtu");
    //VecView(_usr.cg_vars,_viewer);
}

//template <typename REAL>
//void
//ApSolver<REAL>::typeCheck()
//{
//  WxLogger *log = WxLogger::get("apollo-root.console");
//  WxLogStream debStrm = log->getDebugStream();
//  bool allTypeCheck = true;
//  std::ostringstream errorMsg;

//  // loop over each subsolver and make sure all variables being
//  // passed to it are correct
//  typename SubSolverMap_t::iterator ssitr;
//  for (ssitr=_subSolvers.begin(); ssitr!=_subSolvers.end(); ++ssitr)
//  {
//    WxSubSolver<REAL> *ss = ssitr->second;
//    debStrm << "Type checking subsolver '" << ss->name() << "'"
//            << std::endl;

//    // check all read variables
//    std::vector<std::string> readVars = ss->readVarNames();
//    // check if number of variables is correct
//    if (ss->hasVariableReadVars())
//    { // arbitrary number of vars
//      if (ss->numActualReadVars() < ss->numReadVars())
//      {
//        allTypeCheck = false;
//        errorMsg << "Subsolver " << ss->name()
//                 << " expects at least " << ss->numReadVars() << " read variables"
//                 << " but got only " << ss->numActualReadVars()
//                 << std::endl;
//      }
//    }
//    else
//    { // fixed number of vars
//      if (ss->numReadVars() != ss->numActualReadVars())
//      {
//        allTypeCheck = false;
//        errorMsg << "Subsolver " << ss->name()
//                 << " expects " << ss->numReadVars() << " read variables"
//                 << " but got " << ss->numActualReadVars()
//                 << std::endl;

//      }
//    }
//    // check variable types
//    unsigned count = 0;
//    for (unsigned i=0; i<ss->numReadVars(); ++i, ++count)
//    {
//      if (_variables[readVars[i]]->type() != ss->getReadVarType(i))
//      {
//        errorMsg << "Variable " << readVars[i]
//                 << " in subsolver " << ss->name()
//                 << " does not type-check" << std::endl;
//        allTypeCheck = false;
//      }
//    }
//    if (ss->hasVariableReadVars())
//    {
//      for (unsigned i=count; i<ss->numActualReadVars(); ++i)
//      {
//        if (_variables[readVars[i]]->type() != ss->getLastReadVarType())
//        {
//          errorMsg << "Variable " << readVars[i]
//                   << " in subsolver " << ss->name()
//                   << " does not type-check" << std::endl;
//          allTypeCheck = false;
//        }
//      }
//    }

//    // check all write variables
//    std::vector<std::string> writeVars = ss->writeVarNames();
//    // check if number of variables is correct
//    if (ss->hasVariableWriteVars())
//    { // arbitrary number of vars
//      if (ss->numActualWriteVars() < ss->numWriteVars())
//      {
//        allTypeCheck = false;
//        errorMsg << "Subsolver " << ss->name()
//                 << " expects at least " << ss->numWriteVars() << " write variables"
//                 << " but got only " << ss->numActualWriteVars()
//                 << std::endl;
//      }
//    }
//    else
//    { // fixed number of vars
//      if (ss->numWriteVars() != ss->numActualWriteVars())
//      {
//        allTypeCheck = false;
//        errorMsg << "Subsolver " << ss->name()
//                 << " expects " << ss->numWriteVars() << " write variables"
//                 << " but got " << ss->numActualWriteVars()
//                 << std::endl;

//      }
//    }
//    // check variable types
//    count = 0;
//    for (unsigned i=0; i<ss->numWriteVars(); ++i, ++count)
//    {
//      if (_variables[writeVars[i]]->type() != ss->getWriteVarType(i))
//      {
//        errorMsg << "Variable " << writeVars[i]
//                 << " in subsolver " << ss->name()
//                 << " does not type-check" << std::endl;
//        allTypeCheck = false;
//      }
//    }
//    if (ss->hasVariableWriteVars())
//    {
//      for (unsigned i=count; i<ss->numActualWriteVars(); ++i)
//      {
//        if (_variables[writeVars[i]]->type() != ss->getLastWriteVarType())
//        {
//          errorMsg << "Variable " << writeVars[i]
//                   << " in subsolver " << ss->name()
//                   << " does not type-check" << std::endl;
//          allTypeCheck = false;
//        }
//      }
//    }
//  }

//  if (allTypeCheck == false)
//  {
//    WxExcept wxe;
//    wxe << "Type checking error:  " << std::endl
//        << errorMsg.str();
//    throw wxe;
//  }
//}

//template <typename REAL>
//void
//ApSolver<REAL>::writeData(WxIoBase *io, unsigned frame, REAL tcurr, REAL telapsed)
//{

//}

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
        WxStepperStatus<REAL> res = _subSolvers[*ssitr]->step(0.0);
        if (res.getStatus() == false)
        {
          WxExcept wxe("Subsolver ");
          wxe << *ssitr << " failed" << std::endl;
          throw wxe;
        }
      }
    }
}

//template <typename REAL>
//void
//ApSolver<REAL>::writeOnly(WxIoBase *io, WxIoNodeType ioNode)
//{

//}

//template <typename REAL>
//void
//ApSolver<REAL>::endOnly()
//{

//}

template <typename REAL>
void
ApSolver<REAL>::init()
{
    WxLogger *log = WxLogger::get("apollo-root.console");
    WxLogStream debStrm = log->getDebugStream();

    // set current time
    this->setCurrentTime(_tstart);
    // set frame number
    this->setStartFrame(0);

    // initialize subsolvers
    typename SubSolverMap_t::iterator itr;
    for (itr = _subSolvers.begin(); itr != _subSolvers.end(); ++itr)
        itr->second->init();

    // run startOnly subsolvers
    std::string fname = this->runName() + "_0.vtu";
    _outfname = &fname[0];
    PetscViewerFileSetName(_viewer, _outfname);
    debStrm << "Running StartOnly steps...\n" << std::endl;
    startOnly();
}

//template <typename REAL>
//void
//ApSolver<REAL>::load(WxIoBase& io, const WxIoNodeType& solverGrp)
//{
    // Load data to restart ...
  // read in all variables
//  typename VariableMap_t::iterator varItr;
//  for (varItr = _variables.begin(); varItr!=_variables.end(); ++varItr)
//    if ((varItr->second)->getWriteOut())
//      // only load a field if it was written out in the first place
//      (varItr->second)->load(io, solverGrp);

  // load all subsolver data
//  typename SubSolverMap_t::iterator itr;
//  for (itr = _subSolvers.begin(); itr != _subSolvers.end(); ++itr)
//    itr->second->load(io, solverGrp);
//}

template <typename REAL>
void
ApSolver<REAL>::advance(REAL tstart, REAL tend, REAL&dt)
{
    // put ts looping structure here
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
    debStrm << "Partitioning the domain using --> " << _partitioner  << std::endl;
    DMPlexDistribute(*dm, _partitioner, 0, NULL, &dmDist);
    if (dmDist){
        DMDestroy(dm);
        *dm   = dmDist;}
    // get any additional parameters for DM from the command line
    DMSetFromOptions(*dm);
    // SetUp the data structures
    DMSetUp(*dm);
    PetscObjectSetName((PetscObject) *dm, "Mesh");
}

template<typename REAL>
void
ApSolver<REAL>::SetupLocalSpace(DM dm, UserContext usr)
{
    PetscSection   stateSection;
    PetscInt       dof = 1, cStart, cEnd, c;

    DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd);
    PetscSectionCreate(PetscObjectComm((PetscObject)dm), &stateSection);
    PetscSectionSetNumFields(stateSection,1);
    PetscSectionSetFieldComponents(stateSection,0,1);
    PetscSectionSetChart(stateSection, cStart, cEnd);

    for (c = cStart; c < cEnd; ++c)
    {
        PetscSectionSetFieldDof(stateSection,c,0,1);
        PetscSectionSetDof(stateSection, c, dof);
    }

//    for (c = cEndInterior; c < cEnd; ++c)
//        PetscSectionSetConstraintDof(stateSection, c, dof);

//    cind[0] = 0;
//    for (c = cEndInterior; c < cEnd; ++c)
//        PetscSectionSetConstraintIndices(stateSection, c, cind);

    PetscSectionSetUp(stateSection);
    DMSetDefaultSection(dm,stateSection);
    PetscSectionDestroy(&stateSection);
}

// instantiations
template class ApSolver<float>;
template class ApSolver<double>;
