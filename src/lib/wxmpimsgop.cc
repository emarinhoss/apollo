// WarpX lib includes
#include "wxmpimsgop.h"

WxMpiMsgOp::WxMpiMsgOp()
{
  _opMap[WX_MSG_SUM] = MPI_SUM;
  _opMap[WX_MSG_MIN] = MPI_MIN;
  _opMap[WX_MSG_MAX] = MPI_MAX;
  _opMap[WX_MSG_AND] = MPI_LAND;
}

MPI_Op
WxMpiMsgOp::getOp(WxMsgOp op)
{
  return _opMap[op];
}
