#ifndef __wxmpimsgop__
#define __wxmpimsgop__

// WarpX lib includes
#include "wxmsgop.h"

// mpi includes
#include <mpi.h>

// std includes
#include <map>

class WxMpiMsgOp : public WxMsgOpC
{
  public:
    WxMpiMsgOp();

/**
 * Returns an MPI operation corresponding to a string key
 *
 * @param op an operation
 * @return an MPI operation
 */
    MPI_Op getOp(WxMsgOp op);

  private:
    std::map<WxMsgOp, MPI_Op> _opMap;
};

#endif // __wxmpimsgop__
