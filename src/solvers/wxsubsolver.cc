// WarpX lib includes
#include <wxlogger.h>
#include <wxlogstream.h>

// WarpX subsolver includes
#include "wxsubsolver.h"
#include "apsolver.h"

// std includes
#include <sstream>

template <typename REAL>
WxSubSolver<REAL>::WxSubSolver(const std::string& name)
  : WxStepper<REAL>(name), _parent(0), _lastReadVarType(0), _lastWriteVarType(0) 
{
}

template <typename REAL>
WxSubSolver<REAL>::~WxSubSolver()
{
}

template <typename REAL>
const
std::type_info&
WxSubSolver<REAL>::type() const
{ // just return WxSubSolver typeid
  return typeid(WxSubSolver<REAL>);
}

//template <typename REAL>
//WxGridBox<REAL>
//WxSubSolver<REAL>::getGrid() const
//{
//  return this->_parent->getGrid(_onGrid);
//}

template <typename REAL>
void
WxSubSolver<REAL>::setup(const WxCryptSet& wxc)
{
  std::stringstream ss;
  WxLogger *l = WxLogger::get("warpx-root.console");
  WxLogStream debStrm = l->getDebugStream();
  debStrm << "Setting up subsolver '" << wxc.name() << "' of kind '" << this->name() << "'" 
          << std::endl;

  // set names of read and write variables
  if (wxc.has("ReadVars"))
  {
    std::vector<WxAny> readVars = wxc.get<std::vector<WxAny> >("ReadVars");
    std::vector<WxAny>::const_iterator itr;
    for (itr=readVars.begin(); itr!=readVars.end(); ++itr)
      _readVars.push_back( wx_any_cast<std::string>(*itr) );
  }
  if (wxc.has("WriteVars"))
  {
    std::vector<WxAny> writeVars = wxc.get<std::vector<WxAny> >("WriteVars");
    std::vector<WxAny>::const_iterator itr;
    for (itr=writeVars.begin(); itr!=writeVars.end(); ++itr)
      _writeVars.push_back( wx_any_cast<std::string>(*itr) );
  }
  if (this->needsGrid())
  {
    if (wxc.has("OnGrid"))
      _onGrid = wxc.get<std::string>("OnGrid");
    else
    {
      WxExcept wxe;
      wxe << "Subsolver " << this->name()
          << " needs a grid, but OnGrid was not defined in input file";
      throw wxe;
    }
  }
}

template <typename REAL>
void
WxSubSolver<REAL>::declareTypes()
{
}

template <typename REAL>
std::vector<std::string>
WxSubSolver<REAL>::readVarNames() const 
{
  return _readVars;
}

template <typename REAL>
std::vector<std::string>
WxSubSolver<REAL>::writeVarNames() const 
{
  return _writeVars;
}

template <typename REAL>
void
WxSubSolver<REAL>::setReadVarType(unsigned pos, const std::type_info& type)
{
  this->_readVarType[pos] = &type;
}

template <typename REAL>
const std::type_info&
WxSubSolver<REAL>::getReadVarType(unsigned pos)
{
  return *(this->_readVarType[pos]);
}

template <typename REAL>
void
WxSubSolver<REAL>::setWriteVarType(unsigned pos, const std::type_info& type)
{
  this->_writeVarType[pos] = &type;
}

template <typename REAL>
const std::type_info&
WxSubSolver<REAL>::getWriteVarType(unsigned pos)
{
  return *this->_writeVarType[pos];
}

template <typename REAL>
void
WxSubSolver<REAL>::setLastReadVarType(const std::type_info& type)
{
  _lastReadVarType = &type;
}

template <typename REAL>
const std::type_info&
WxSubSolver<REAL>::getLastReadVarType()
{
  return *this->_lastReadVarType;
}

template <typename REAL>
void
WxSubSolver<REAL>::setLastWriteVarType(const std::type_info& type)
{
  _lastWriteVarType = &type;
}

template <typename REAL>
const std::type_info&
WxSubSolver<REAL>::getLastWriteVarType()
{
  return *this->_lastWriteVarType;
}

template <typename REAL>
unsigned
WxSubSolver<REAL>::numReadVars() const
{
  return this->_readVarType.size();
}

template <typename REAL>
unsigned
WxSubSolver<REAL>::numWriteVars() const
{
  return this->_writeVarType.size();
}

template <typename REAL>
bool
WxSubSolver<REAL>::hasVariableReadVars() const
{
  return this->_lastReadVarType ? true : false;
}

template <typename REAL>
bool
WxSubSolver<REAL>::hasVariableWriteVars() const
{
  return this->_lastWriteVarType ? true : false;
}

template <typename REAL>
bool
WxSubSolver<REAL>::needsGrid() const
{
  return true;
}

template <typename REAL>
void
WxSubSolver<REAL>::setParent(ApSolver<REAL> *parent)
{
  _parent = parent;
}

template <typename REAL>
ApSolver<REAL> *WxSubSolver<REAL>::getParent() const
{
  return _parent;
}

template <typename REAL>
unsigned
WxSubSolver<REAL>::numActualReadVars() const
{
  return this->_readVars.size();
}

template <typename REAL>
unsigned
WxSubSolver<REAL>::numActualWriteVars() const 
{
  return this->_writeVars.size();
}

template <typename REAL>
std::vector<std::string>
WxSubSolver<REAL>::getWriteVarList() const 
{
  return _writeVars;
}

template <typename REAL>
std::vector<std::string>
WxSubSolver<REAL>::getReadVarList() const
{
  return _readVars;
}

template <typename REAL>
void
WxSubSolver<REAL>::setWriteVarList(const std::vector<std::string>& wr) 
{
  _writeVars.erase(_writeVars.begin(), _writeVars.end());
  std::vector<std::string>::const_iterator itr;
  for (itr=wr.begin(); itr!=wr.end(); ++itr)
    _writeVars.push_back( *itr );
}

template <typename REAL>
void
WxSubSolver<REAL>::setReadVarList(const std::vector<std::string>& rv)
{
  _readVars.erase(_readVars.begin(), _readVars.end());
  std::vector<std::string>::const_iterator itr;
  for (itr=rv.begin(); itr!=rv.end(); ++itr)
    _readVars.push_back( *itr );
}

// instantiations
template class WxSubSolver<float>;
template class WxSubSolver<double>;
