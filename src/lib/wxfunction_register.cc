// WarpX includes
#include "wxcreator.h"

#include "wxfunction.h"
#include "wxconstfunc.h"
#include "wxconstvecfunc.h"
#include "wxcombofunc.h"
#include "wxexprfunction.h"

// register constant function
WxCreator< WxConstFunc<float>, WxFunction<float> > __constFunc_f("constFunc");
WxCreator< WxConstFunc<double>, WxFunction<double> > __constFunc_d("constFunc");

WxCreator< WxConstFunc<int, float>, WxFunction<int, float> > __constFunc_if("constFunc");
WxCreator< WxConstFunc<int, double>, WxFunction<int, double> > __constFunc_id("constFunc");

// register constant vector function
WxCreator< WxConstVecFunc<float>, WxFunction<float> > __constVecFunc_f("constVecFunc");
WxCreator< WxConstVecFunc<double>, WxFunction<double> > __consVecFunc_d("constVecFunc");

WxCreator< WxConstVecFunc<int, float>, WxFunction<int, float> > __constVecFunc_if("constVecFunc");
WxCreator< WxConstVecFunc<int, double>, WxFunction<int, double> > __consVecFunc_id("constVecFunc");

// register combo-function function
WxCreator< WxComboFunc<float>, WxFunction<float> > __comboFunc_f("comboFunc");
WxCreator< WxComboFunc<double>, WxFunction<double> > __comboFunc_d("comboFunc");

// register expression evaluator function
WxCreator< WxExprFunc, WxFunction<double> > __exprfunc_d("exprFunc");
WxCreator< WxExprFunc, WxFunction<float> > __exprfunc_f("exprFunc");
