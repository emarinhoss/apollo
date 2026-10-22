// WarpX lib includes
#include <wxexprparser.h>
#include <wxexcept.h>

// std includes
#include <algorithm>
#include <set>

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

  checkForUndefinedNames();
}

void
WxExprParser::checkForUndefinedNames()
{
  // The variable factory above hands muParser a fresh zero for any identifier
  // it has never seen. That is what lets a pre-expression introduce an
  // intermediate ("r = sqrt(...)") without declaring it first - but it is also
  // why a misspelt constant evaluates to 0.0 with no diagnostic at all. On the
  // shipped isentropic vortex deck, changing one character of "beta" to "bta"
  // in the velocity expression shifts the initial x-momentum by 47% and the
  // energy by 29%, and the run still exits 0.
  //
  // muParser cannot distinguish "assigning a new intermediate" from "reading a
  // name that does not exist", so the check happens here instead: every
  // identifier the expression actually uses must be an independent variable, a
  // declared constant, one of our result temporaries, or assigned by one of the
  // pre-expressions. Anything else was invented by the factory.

  // Names introduced by assignment in a pre-expression, e.g. "u = uo - ...".
  std::set<std::string> assigned;
  for (unsigned i=0; i<preExprs.size(); ++i)
  {
    const std::string& e = preExprs[i];
    for (std::string::size_type k=0; k<e.size(); ++k)
    {
      if (e[k] != '=') continue;
      // Skip the comparison operators: ==, <=, >=, !=
      if (k+1 < e.size() && e[k+1] == '=') break;
      if (k > 0 && (e[k-1]=='<' || e[k-1]=='>' || e[k-1]=='!' || e[k-1]=='=')) break;
      std::string lhs = e.substr(0, k);
      // trim
      std::string::size_type b = lhs.find_first_not_of(" \t");
      std::string::size_type d = lhs.find_last_not_of(" \t");
      if (b != std::string::npos)
        assigned.insert(lhs.substr(b, d-b+1));
      break; // only the first '=' is the assignment
    }
  }

  const mu::varmap_type& used = parser.GetUsedVar();
  std::vector<std::string> unknown;
  for (mu::varmap_type::const_iterator i=used.begin(); i!=used.end(); ++i)
  {
    const std::string& name = i->first;
    if (assigned.count(name))                       continue;
    if (consts.find(name) != consts.end())          continue;
    if (std::find(varNames.begin(), varNames.end(), name) != varNames.end()) continue;
    // Our own result temporaries, __res0__ and friends.
    if (name.size() > 4 && name.compare(0, 5, "__res") == 0) continue;
    unknown.push_back(name);
  }

  if (!unknown.empty())
  {
    WxExcept wxe("Undefined name");
    wxe << (unknown.size() > 1 ? "s " : " ");
    for (unsigned i=0; i<unknown.size(); ++i)
      wxe << (i ? ", " : "") << "'" << unknown[i] << "'";
    wxe << " in expression. Every name must be a coordinate (";
    for (unsigned i=0; i<varNames.size(); ++i)
      wxe << (i ? ", " : "") << varNames[i];
    wxe << "), a constant declared in this block, or assigned by an earlier "
        << "entry of progn. Check the spelling, and note that a constant "
        << "written without a decimal point is read as an integer and is not "
        << "exported to expressions." << std::endl;
    throw wxe;
  }
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
