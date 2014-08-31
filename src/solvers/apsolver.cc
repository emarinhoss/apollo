// WarpX lib includes
#include <wxcreator.h>
#include <wxlogger.h>
#include <wxlogstream.h>
#include <wxtimer.h>

// WarpX solver includes
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
    WxLogStream debStrm = WxLogger::get("warpx-root.console")->getDebugStream();

    // setup our parent first
    WxSolverBase<REAL>::setup(wxc);
    // start and end time for simulation
    std::vector<WxAny> times = wxc.template get<std::vector<WxAny> >("Time");
    _tstart = wx_any_cast<REAL>(times[0]);
    _tend = wx_any_cast<REAL>(times[1]);

    // number of files to write
    _nout = wxc.template get<int>("Out");

    // time step to take
    _dt = wxc.template get<REAL>("Dt");
}

template <typename REAL>
WxSubSolver<REAL>*
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

  REAL temp_var;

  unsigned nout;
  REAL tsize; // time between file output

  // write data to file before running main loop
  //this->writeData(&this->getIo(), frame, tcurr, 0.0);

  // main solver loop

  frame += 1;
  nout = _nout; //not sure how to take restart into account yet

  REAL timeStep = _dt; // initial time step to use
  for (unsigned i=0; i<nout; ++i)
  {
//    infStrm << "Advancing solution"
//            << " from time " << tcurr
//            << " to " << temp_var
//            << "..."
//            << std::endl;

    WxTimer advTimer;
    // advance solution on each block by 'tsize'
    advTimer.startTimer();
    advance(tcurr, temp_var, timeStep);
    advTimer.stopTimer();

    // write solution to file
//    this->writeData(&this->getIo(), frame, temp_var, advTimer.secondsElapsed());

//    infStrm << "Advance completed in "
//            << advTimer.timeElapsedAsString()
//            << std::endl << std::endl;

    // advance tcurr and frame number
//    tcurr += tsize;
//    frame += 1;
  }
}

template <typename REAL>
WxStepperStatus<REAL>
ApSolver<REAL>::step(REAL dt)
{
  return WxStepperStatus<REAL>();
}

//template <typename REAL>
//void
//ApSolver<REAL>::typeCheck()
//{
//  WxLogger *log = WxLogger::get("warpx-root.console");
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

//template <typename REAL>
//void
//ApSolver<REAL>::startOnly()
//{

//}

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
    // Initialize subsolvers

  // set current time
  this->setCurrentTime(_tstart);
  // set frame number
  this->setStartFrame(0);

  // initialize subsolvers
  typename SubSolverMap_t::iterator itr;
  for (itr = _subSolvers.begin(); itr != _subSolvers.end(); ++itr)
    itr->second->init();

  // run startOnly subsolvers
  // debStrm << "Running StartOnly steps...\n" << std::endl;
  // startOnly();
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

// instantiations
template class ApSolver<float>;
template class ApSolver<double>;
