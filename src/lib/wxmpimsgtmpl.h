#ifndef __wxmpimsgtmpl__
#define __wxmpimsgtmpl__

// WarpX includes
#include "wxmsgtmpl.h"
#include "wxmpitraits.h"
#include "wxmpimsgop.h"

// std includes
#include <iostream>

/**
 * Mpi specific message status wrapper
 */
struct WxMpiMsgStatus_v : public WxMsgStatus_v
{
    WxMpiMsgStatus_v()
            : data(0) {
    }
    MPI_Request request; // send/recv request
    void *data; // pointer to memory for recv
};

template<typename T>
class WxMpiMsgTmpl : public WxMsgTmpl<T>
{
  public:
    WxMpiMsgTmpl(MPI_Comm comm)
            : _comm(comm), _sending(false) {
    }

/**
 * Send an array to another rank.
 *
 * @param arr std::vector of data being sent
 * @param recvRank rank that will receive array
 */
    void send(const std::vector<T>& arr, unsigned recvRank, int tag);

/**
 * Send an array to another rank.
 * 
 * @parem num number of elements to send
 * @param arr array of length 'num' of data being sent
 * @param recvRank rank that will receive array
 */
    void send(unsigned num, T* arr, unsigned recvRank, int tag);

/**
 * Receive an array from another rank.
 *
 * @param num number of elements to reciev
 * @param array array that is filled with received values
 * @param sendRank rank that sent array
 */
    void recv(unsigned num, std::vector<T>& array, unsigned sendRank, int tag);

/**
 * Receive an array from another rank.
 *
 * @param num number of elements to receive
 * @param array array that is filled with received values
 * @param sendRank rank that sent array
 */
    void recv(unsigned num, T* array, unsigned sendRank, int tag);

/**
 * Receive an array from another rank. This is a non-blocking call.
 * 
 * @param num number of elements to receive
 * @param array array that is filled with received values
 * @param sendRank rank that sent array
 * @return message status
 */
    WxMsgStatus startRecv(unsigned num, unsigned sendRank, int tag);

/**
 * Reduce data to all ranks
 * 
 * @param num Number of elements being reduced
 * @param sendBuff buffer to reduce
 * @param recvBuff buffer to recieve reduced data
 * @param op operation to perform. These are one of those listed in WxMsgOp enum
 */
    void allReduce(unsigned num, T* sendBuff, T* recvBuff, WxMsgOp op);

  private:
    MPI_Comm _comm;
    bool _sending;
    MPI_Request _sendRequest;
    WxMpiMsgOp _ops;
};

#endif //  __wxmpimsgtmpl__
