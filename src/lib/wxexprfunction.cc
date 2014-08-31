// WarpX lib includes
#include "wxexprfunction.h"

void
WxExprFunc::setup(const WxCryptSet& wxc) 
{
  // set t, x, y, z as independent variables
  _parser.appendIndVar("t");
  _parser.appendIndVar("x");
  _parser.appendIndVar("y");
  _parser.appendIndVar("z");
  // read in constants
  std::vector<std::string> dbls = wxc.getNames<double>();
  for (unsigned i=0; i<dbls.size(); ++i)
  {
    double val = wxc.get<double>(dbls[i]);
    _parser.addConstant(dbls[i], val);
  }
  // read expressions in progn list
  if (wxc.has("progn")) 
  {
    std::vector<std::string> progn 
      = wxc.getVec<std::string>("progn");
    for (unsigned i=0; i<progn.size(); ++i)
      _parser.appendPreExpr(progn[i]);
  }
  // read exprList
  std::vector<std::string> exprList 
    = wxc.getVec<std::string>("exprList");
  for (unsigned i=0; i<exprList.size(); ++i)
  {
    unsigned ec = _parser.addExpr(exprList[i]);
    // put result code into a vector
    _resCodes.push_back(ec);
  }
 // these many components will be set
  _nres = exprList.size();
  // setup parser
  _parser.setup();
}

bool
WxExprFunc::func(unsigned n, double *tx, double *r) 
{
  std::vector<double> inp(4);
  // set variables for use in expressions
  inp[0] = tx[0];
  if (n>1)
    inp[1] = tx[1];
  if (n>2)
    inp[2] = tx[2];
  if (n>3)
    inp[3] = tx[3];
  _parser.eval(inp);
  // set all the result 
    for (unsigned j=0; j<_nres; ++j)
      r[j] = _parser.result(_resCodes[j]);

  return true;
}
