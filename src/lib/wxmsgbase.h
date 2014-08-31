#ifndef __wxmsgbase__
#define __wxmsgbase__

// WarpX includes
#include "wxmsgtmpl.h"
#include "wxtypelist.h"
#include "wxdatatypes.h"
#include "wxexcept.h"

// std includes
#include <vector>

/**
 * Provides an abstract interface for message based communication
 * between different processes.
 */
class WxMsgBase 
{
  public:

/**
 * Destructor
 */
    virtual ~WxMsgBase() {
    }

/**
 * Rank of process
 * 
 * @return this rank
 */
    virtual int rank() const = 0;

/**
 * Number of processes taking part in messaging
 *
 * @return num of processes
 */
    virtual unsigned numProcs() const = 0;

/**
 * Get parent communicating processor group
 *
 * @return parent proccesor group
 */
    WxMsgBase* parent() const {
      return _parent;
    }

/**
 * Split communicator into a child communicator
 *
 * @param ranks list of processors in old communicator which are to be
 * in the new communicator
 * @return new communicator
 */
    virtual WxMsgBase* createSubComm(const std::vector<int>& ranks) = 0;


/**
 * Block till all processes hit this barrier
 */
    virtual void barrier() const = 0;

/**
 * Send a std::vector to another rank.
 * 
 * @param array std::vector of data being sent
 * @param recvRank rank that will receive array
 */
    template<typename T>
    void send(const std::vector<T>& array, unsigned recvRank, int tag=-1) {
      this->_getMsg<T>()->send(array, recvRank, (tag == -1) ? _sendTag : tag);
    }

/**
 * Send an array to another rank.
 * 
 * @parem num number of elements to send
 * @param array array of length 'num' of data being sent
 * @param recvRank rank that will receive array
 */
    template<typename T>
    void send(unsigned num, T* array, unsigned recvRank, int tag=-1) {
      this->_getMsg<T>()->send(num, array, recvRank, (tag == -1) ? _sendTag : tag);
    }

/**
 * Receive a std::vector from another rank.
 * 
 * @param num number of elements to receive
 * @param array array that is filled with received values
 * @param sendRank rank that sent array
 */
    template<typename T>
    void recv(int num, std::vector<T>& array, unsigned sendRank, int tag=-1) {
      this->_getMsg<T>()->recv(num, array, sendRank, (tag == -1) ? _recvTag : tag);
    }

/**
 * Receive an array from another rank.
 * 
 * @param num number of elements to receive
 * @param array array that is filled with received values
 * @param sendRank rank that sent array
 */
    template<typename T>
    void recv(unsigned num, T* array, unsigned sendRank, int tag=-1) {
      this->_getMsg<T>()->recv(num, array, sendRank, (tag == -1) ? _recvTag : tag);
    }

/**
 * Receive an array from another rank. This is a non-blocking call.
 * 
 * @param num number of elements to receive
 * @param array array that is filled with received value
 * @param sendRank rank we want to receive from
 * @return message status
 */
    template<typename T>
    WxMsgStatus startRecv(unsigned num,  unsigned sendRank, int tag=-1) {
      return this->_getMsg<T>()->startRecv(
          num, sendRank, (tag == -1) ? _recvTag : tag);
    }

/**
 * Finish the receive started by startRecv and return a pointer to the
 * data recieved. The calling function owns the pointer and is
 * resposible for freeing it.
 *
 * @param ms message status object returned by startRecv
 * @return pointer to the data recieved. Caller owns the returned pointer
 */
    virtual void * finishRecv(WxMsgStatus ms) = 0;

/**
 * Check status of recieve started by a startRecv. If this call
 * returns true then finishRecv can be called to immediately recieve
 * the data.
 *
 * @param ms message status object returned by startRecv
 * @return true it recieve has been completed, false otherwise
 */
    virtual bool checkRecv(WxMsgStatus ms) = 0;

/**
 * Reduce data to all ranks
 * 
 * @param num Number of elements being reduced
 * @param sendBuff buffer to reduce
 * @param recvBuff buffer to recieve reduced data
 * @param op operation to perform. These are one of those listed in WxMsgOp enum
 */
    template<typename T>
    void allReduce(unsigned num, T* sendBuff, T* recvBuff, WxMsgOp op) {
      this->_getMsg<T>()->allReduce(num, sendBuff, recvBuff, op);
    }

  protected:

/**
 * Constructor. This is protected so only children can make instances.
 */
    WxMsgBase(int sendTag=0, int recvTag=0, WxMsgBase* parent=0)
            : _sendTag(sendTag), _recvTag(recvTag), _parent(parent) {
    }

/**
 * Add a new messager : the derived class should call this to setup
 * WxMsgBase properly
 */
    template <typename T>
    void addMsg(WxMsgTmpl<T> *b) {
      wxTypeMapExtract<T>(_msgTypeMap)._msg = b;
    }

  private:

    // To prevent use
    WxMsgBase(const WxMsgBase&);
    WxMsgBase& operator=(const WxMsgBase&);

    int _sendTag, _recvTag;
    WxMsgBase *_parent;

/**
 * Get a messager object with the proper type
 */
    template <typename T>
    WxMsgTmpl<T>* _getMsg() {
      WxMsgTmpl<T> *r = wxTypeMapExtract<T>(_msgTypeMap)._msg;
      if (r) return r;
      WxExcept wxe;
      wxe << "Message type not set properly"; 
      throw wxe;
    }

  public:
    //  xlC doesn't like internal private structs

    // Container class for all message-ers
    template<typename T>
    struct WxMsgContainer {
        WxMsgContainer() : _msg(0) {
        }
        virtual ~WxMsgContainer() {
          delete _msg;
        }
        // this points to a derived class of WxMsgTmpl<T>
        WxMsgTmpl<T> *_msg;
    };

    // Objects of type WxMsgTypeMap_t inherit from all WxMsgTmpl<T>
    // where T belongs to the WxMsgTypelist_t. Thus it acts like a
    // container for all message-er objects in the system.
    typedef WxTypeMap<WxDataTypes_t, WxMsgContainer> WxMsgTypeMap_t;

    WxMsgTypeMap_t _msgTypeMap; // container of communicators
};

#endif //  __wxmsgbase__
