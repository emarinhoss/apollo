#ifndef __wxselfmsgtmpl__
#define __wxselfmsgtmpl__

// WarpX includes
#include "wxmsgtmpl.h"

// std includes
#include <iostream>

/**
 * Self communicator specific message status wrapper
 */
struct WxSelfMsgStatus_v : public WxMsgStatus_v
{
};

template<typename T>
class WxSelfMsgTmpl : public WxMsgTmpl<T>
{
  public:
    WxSelfMsgTmpl() {
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
};

#endif //  __wxselfmsgtmpl__

