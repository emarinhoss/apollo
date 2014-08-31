// WarpX lib includes
#include <wxexprparser.h>
#include <wxexcept.h>

mu::value_type
Begin(const mu::value_type* a_afArg, int a_iArgc) 
{
  mu::value_type fRes=0;
  return fRes;
}

double*
addVariable(const mu::char_type* a_szName, void* a_pUserData) 
{
  // at most MAX_VARS variables can be defined
#define MAX_VARS 5000
  static double afValBuf[MAX_VARS];
  static int iVal = 0;

  afValBuf[iVal] = 0;
  if (iVal>=(MAX_VARS-1))
    throw mu::ParserError( _T("Variable buffer overflow.") );
  return& afValBuf[iVal++];
}

WxExprParser::WxExprParser()
  : resCount(0) 
{
  // set variable construction factory
  parser.SetVarFactory(addVariable, &parser);
  // set expression sequence function
  parser.DefineFun(_T("begin"), Begin, false);
}

void
WxExprParser::appendIndVar(const std::string& v) 
{
  varNames.push_back(v);
  vars.push_back(0.0);
}

void
WxExprParser::addConstant(const std::string& name, double value) 
{
  consts[name] = value;
}

void
WxExprParser::appendPreExpr(const std::string& expr) 
{
  preExprs.push_back(expr);
}

unsigned
WxExprParser::addExpr(const std::string& expr) 
{
  unsigned count = resCount;
  exprs.push_back(expr);  // add expression
  res.push_back(0.0); // set its result value to 0.0
  resCount++;
  return count;
}

void
WxExprParser::setup() 
{
  std::ostringstream exprStrm;

  // set independent variables
  for (unsigned i=0; i<vars.size(); ++i)
    parser.DefineVar(varNames[i], &vars[i]);

  // add constants
  std::map<std::string, double>::iterator citr;
  for (citr=consts.begin(); citr!=consts.end(); ++citr)
    parser.DefineVar(citr->first, &citr->second);

  // Construct string representing expression. This is of the form
  //
  // begin(preExpr, preExpr, .., __res0__ = expr, __res1__ = expr, ..)
  //
  // The begin() form ensures that the evaluation occurs from left to
  // right.
  exprStrm << "begin(";
  // add pre-expressions in order they were appended
  unsigned si;
  for (si=0; si<preExprs.size(); ++si) {
    exprStrm << preExprs[si] << ", ";
  }
  // add expressions
  for (si=0; si<exprs.size()-1; ++si) {
    std::ostringstream vn;
    vn << "__res" << si << "__";  // temp. name for result
    parser.DefineVar(vn.str(), &res[si]);
    exprStrm << vn.str() << "=" << exprs[si] << ", ";
  }
  std::ostringstream vn;
  vn << "__res" << si << "__";
  parser.DefineVar(vn.str(), &res[si]);
  exprStrm << vn.str() << "=" << exprs[si] << ")";

  // set expression to parse
  parser.SetExpr(exprStrm.str());
  // copy this string for use in debugging
  fullExprString = exprStrm.str();
}

void
WxExprParser::eval(const std::vector<double>& iv) 
{
  // set independent variables
  for (unsigned i=0; i<iv.size(); ++i)
    vars[i] = iv[i];
  // evaluate expression
  try
  {
    parser.Eval();
  }
  catch (mu::ParserError& e)
  {
    WxExcept wxe("WxExprParser::eval : Error in expression '");
    wxe << e.GetExpr() << "' with message '" << e.GetMsg() << "'" << std::endl;
    throw wxe;
  }
}

double
WxExprParser::result(unsigned rs) const 
{
  return res[rs];
}

std::string
WxExprParser::getExprString() const 
{
  return fullExprString;
}
