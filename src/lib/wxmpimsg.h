#ifndef __wxmpimsg__
#define __wxmpimsg__

#ifdef _DO_USE_MPI_
# include <mpi.h>
#endif

// WarpX includes
#include "wxmsgbase.h"
#include "wxmpimsgtmpl.h"
#include "wxdatatypes.h"

// std includes

class WxMpiMsg : public WxMsgBase
{
  public:

/**
 * Construct a new MPI messenger give a set of communicating
 * processors.
 *
 * @param commProcs set of communicating processors
 */
    WxMpiMsg();

/**
 * Rank of process
 * 
 * @return this rank
 */
    int rank() const 
    {
      int r;
      MPI_Comm_rank(_comm, &r);
      return r;
    }

/**
 * Number of processes taking part in messaging
 *
 * @return num of processes
 */
    unsigned numProcs() const 
    {
      int np;
      MPI_Comm_size(_comm, &np);
      return np;
    }

/**
 * Split communicator into a child communicator
 *
 * @param ranks list of processors in old communicator which are to be
 * in the new communicator
 * @return new communicator
 */
    WxMsgBase* createSubComm(const std::vector<int>& ranks);


/**
 * Block till all processes hit this barrier
 */
    void barrier() const 
    {
      MPI_Barrier(_comm);
    }

/**
 * Finish the receive started by startRecv and return a pointer to the
 * data recieved. The calling function owns the pointer and is
 * resposible for freeing it.
 *
 * @param ms message status object returned by startRecv
 */
    void * finishRecv(WxMsgStatus ms);

/**
 * Check status of recieve started by a startRecv. If this call
 * returns true then finishRecv can be called to immediately recieve
 * the data.
 *
 * @param ms message status object returned by startRecv
 * @return true it recieve has been completed, false otherwise
 */
    bool checkRecv(WxMsgStatus ms);

/**
 * Get MPI communicator for this class
 *
 * @return MPI communicator
 */
    MPI_Comm getMpiComm() 
    {
      return _comm;
    }

  private:
/**
 * Private ctor: for internal use only
 */
    WxMpiMsg(WxMpiMsg *parent, MPI_Comm comm);

    MPI_Comm _comm;
};

#endif //  __wxmpimsg__
