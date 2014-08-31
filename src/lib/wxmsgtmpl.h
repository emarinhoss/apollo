#ifndef __wxmsgtmpl__
#define __wxmsgtmpl__

// WarpX lib includes
#include "wxmsgop.h"

// std includes
#include <vector>

/**
 * Provides a means for derived messengers to return implimentation
 * specific message status flags and data. The returned type is opaque
 * and should not be directly fiddled around with.
 */
struct WxMsgStatus_v {};
typedef WxMsgStatus_v* WxMsgStatus;

template<typename T>
class WxMsgTmpl {
  public:
/**
 * Destructor
 */
    virtual ~WxMsgTmpl() {
        delete [] _sendBuff;
        delete [] _recvBuff;
    }

/**
 * Send an array to another rank.
 *
 * @param array std::vector of data being sent
 * @param recvRank rank that will receive array
 */
    virtual void send(const std::vector<T>& array, unsigned recvRank, int tag) = 0;

/**
 * Send an array to another rank.
 * 
 * @parem num number of elements to send
 * @param array array of length 'num' of data being sent
 * @param recvRank rank that will receive array
 */
    virtual void send(unsigned num, T* array, unsigned recvRank, int tag) = 0;

/**
 * Receive an array from another rank.
 * 
 * @param num number of elements to recieve
 * @param array array that is filled with received values
 * @param sendRank rank that sent array
 */
    virtual void recv(unsigned num, std::vector<T>& array, unsigned sendRank, int tag) = 0;

/**
 * Receive an array from another rank.
 * 
 * @param num number of elements to receive
 * @param array array that is filled with received values
 * @param sendRank rank that sent array
 */
    virtual void recv(unsigned num, T* array, unsigned sendRank, int tag) = 0;

/**
 * Receive an array from another rank. This is a non-blocking call.
 * 
 * @param num number of elements to receive
 * @param sendRank rank that sent array
 * @return message status
 */
    virtual WxMsgStatus startRecv(unsigned num, unsigned sendRank, int tag=-1) = 0;

/**
 * Reduce data to all ranks
 * 
 * @param num Number of elements being reduced
 * @param sendBuff buffer to reduce
 * @param recvBuff buffer to recieve reduced data
 * @param op operation to perform. These are one of those listed in WxMsgOp enum
 */
    virtual void allReduce(unsigned num, T* sendBuff, T* recvBuff, WxMsgOp op) = 0;

  protected:

/**
 * Protected so only children can make instances
 */
    WxMsgTmpl() 
            : _sendSize(0), _sendBuff(0), _recvSize(0), _recvBuff(0) {
    }

/**
 * Resizes a array.  Sets number of elements to the exact number of
 * elements in the data.
 *
 * @param array the array to be resized.
 * @param numElem the number of elements in the array. Is set to reqNumElem
 * @param reqNumElem the number of array elements needed
 */
    void resizeArray(T* &array, unsigned &numElem, unsigned reqNumElem) const {
        if (numElem < reqNumElem) {
            // delete old array and allocate new memory
            delete[] array;
            array = new T[reqNumElem];
        }
        // set number of elements
        numElem = reqNumElem;
        return;
    }

/**
 * Checks the buffer and resizes it if necessary
 * 
 * @param reqSize the needed size for the buffer
 */
    void resizeSendBuff(unsigned reqSize) {
        if (_sendSize < reqSize) {
            // delete old array and allocate new memory
            delete [] _sendBuff;
            _sendSize = reqSize;
            _sendBuff = new T[_sendSize];
        }
        return;
    }

/**
 * Checks the buffer and resizes it if necessary
 *
 * @param reqSize the needed size for the buffer
 */
    void resizeRecvBuff(unsigned reqSize) {
        if (_recvSize < reqSize) {
            delete [] _recvBuff;
            _recvSize = reqSize;
            _recvBuff = new T[_recvSize];
        }
        return;
    }

    unsigned _sendSize;
    T* _sendBuff; // buffer for sending

    unsigned _recvSize;
    T* _recvBuff; // buffer for receiving
};

#endif // __wxmsgtmpl__
