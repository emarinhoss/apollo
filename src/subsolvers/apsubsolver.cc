// WarpX lib includes
#include <wxlogger.h>
#include <wxlogstream.h>

// WarpX subsolver includes
#include "apsubsolver.h"
#include "apsolver.h"

// std includes
#include <sstream>

template <typename REAL>
ApSubSolver<REAL>::ApSubSolver(const std::string& name)
  : WxStepper<REAL>(name), _parent(0)
{
}

template <typename REAL>
ApSubSolver<REAL>::~ApSubSolver()
{
}

template <typename REAL>
const
std::type_info&
ApSubSolver<REAL>::type() const
{ // just return WxSubSolver typeid
  return typeid(ApSubSolver<REAL>);
}

template <typename REAL>
void
ApSubSolver<REAL>::setup(const WxCryptSet& wxc, DM dm)
{
//  std::stringstream ss;
  WxLogger *l = WxLogger::get("apollo-root.console");
  WxLogStream debStrm = l->getDebugStream();
  debStrm << "Setting up subsolver '" << wxc.name() << "' of kind '" << this->name() << "'" 
          << std::endl;

}

template <typename REAL>
void
ApSubSolver<REAL>::setParent(ApSolver<REAL> *parent)
{
  _parent = parent;
}

template <typename REAL>
ApSolver<REAL> *ApSubSolver<REAL>::getParent() const
{
  return _parent;
}

//template <typename REAL>
//std::vector<WxAny> ApSubSolver<REAL>::getDataStructure()
//{
//}

// instantiations
template class ApSubSolver<float>;
template class ApSubSolver<double>;
