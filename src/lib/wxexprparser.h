#ifndef __wxexprparser__
#define __wxexprparser__

// wxetc includes
#include <mup/muParser.h>

// std includes
#include <iostream>
#include <string>
#include <map>
#include <set>
#include <vector>
#include <sstream>

/** Variable factory */
double* addVariable(const mu::char_type* a_szName, void* a_pUserData);
/** Sequence of expressions */
mu::value_type Begin(const mu::value_type* a_afArg, int a_iArgc);

/**
 * A class implementing interfaces to the muParser package. This class
 * is not templated as muParser works with doubles only
 */
class WxExprParser {
  public:

/** Default ctor */
    WxExprParser();

/**
 * Append an indepent variable to system. These are expected to be in
 * the same sequence as given to the eval() function
 */
    void appendIndVar(const std::string& v);

/**
 * Add a new constant
 *
 * @param name Name of constant
 * @param value Value of constant
 */
    void addConstant(const std::string& name, double value);

/**
 * Append a pre-expression: these are evaluated in the sequence in
 * which they are added
 *
 * @param expr Expression string
 */
    void appendPreExpr(const std::string& expr);

/**
 * Add an expression. This returns an integer which can be used to
 * query the results once the system is evaluated.
 *
 * @param expr Expression string
 * @return integer representing expression
 */
    unsigned addExpr(const std::string& expr);

/**
 * Setup the system once constants, pre-expressions and expressions
 * have been added.
 */
    void setup();

/**
 * Evaluate the system of expressions. Results can be queries using
 * the result() function.
 *
 * @param iv Values of indepenent variables in system
 */
    void eval(const std::vector<double>& iv);

/**
 * Return value of expression after evaluation
 *
 * @param rs Expression number returned by addExpr() method
 * @return result of evaluation
 */
    double result(unsigned rs) const;

/**
 * Return string representation of the expression. This is only useful
 * for debugging and generally does not make too much sense.
 *
 * @return string representing expression
 */
    std::string getExprString() const;

  private:
/**
 * Reject identifiers that the variable factory invented. Called at the end of
 * setup(), so a misspelt name is an error before the run starts rather than a
 * silent zero in the middle of it.
 */
    void checkForUndefinedNames();

/** Counter for results */
    unsigned resCount;
/** Parser object */
    mu::Parser parser;
/** Names of indpendent variables */
    std::vector<std::string> varNames;
/** List of independent variables */
    std::vector<double> vars;
/** Constants */
    std::map<std::string, double> consts;
/** Pre-expressions */
    std::vector<std::string> preExprs;
/** Expression to execute */
    std::vector<std::string> exprs;
/** Results map */
    std::vector<double> res;
/** Fully constructed expression (for debugging only) */
    std::string fullExprString;
};

#endif // __wxexprparser__
