#ifndef __wxselfmsg__
#define __wxselfmsg__

// WarpX includes
#include "wxmsgbase.h"
#include "wxselfmsgtmpl.h"
#include "wxdatatypes.h"

// std includes

class WxSelfMsg : public WxMsgBase
{
  public:

/**
 * Construct a new MPI messenger give a set of communicating
 * processors.
 *
 * @param commProcs set of communicating processors
 */
    WxSelfMsg();

/**
 * Rank of process
 * 
 * @return this rank
 */
    int rank() const {
      return 0;
    }

/**
 * Number of processes taking part in messaging
 *
 * @return num of processes
 */
    unsigned numProcs() const {
      return 1;
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
    void barrier() const {
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
};

#endif //  __wxmpimsg__
