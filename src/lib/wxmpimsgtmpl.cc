// WarpX lib includes
#include "wxtypelist.h"
#include "wxdatatypes.h"
#include "wxmpimsgtmpl.h"

template<typename T>
void 
WxMpiMsgTmpl<T>::send(const std::vector<T>& arr, unsigned recvRank, int tag) 
{
  unsigned num = arr.size();
  this->resizeSendBuff(num);
  for (unsigned i=0; i<num; ++i)
    this->_sendBuff[i] = arr[i];

  MPI_Send( (void*) this->_sendBuff,
            num,
            WxMpiTraits<T>::mpiType(),
            recvRank,
            tag,
            _comm );
}

template<typename T>
void 
WxMpiMsgTmpl<T>::send(unsigned num, T* arr, unsigned recvRank, int tag) 
{
  MPI_Send( (void*) arr,
            num,
            WxMpiTraits<T>::mpiType(),
            recvRank,
            tag,
            _comm );
}

template<typename T>
void 
WxMpiMsgTmpl<T>::recv(unsigned num, std::vector<T>& array, unsigned sendRank, int tag) 
{
  MPI_Status status;
  this->resizeRecvBuff(num);
  // do a blocking receive
  MPI_Recv( (void*) this->_recvBuff,
            num,
            WxMpiTraits<T>::mpiType(),
            sendRank,
            tag,
            _comm,
            &status );
  // stick stuff into the vector
  for (unsigned i=0; i<num; ++i)
    array.push_back(this->_recvBuff[i]);
}

template<typename T>
void 
WxMpiMsgTmpl<T>::recv(unsigned num, T* array, unsigned sendRank, int tag) 
{
  MPI_Status status;
  // do a blocking receive
  MPI_Recv( (void*) array,
            num,
            WxMpiTraits<T>::mpiType(),
            sendRank,
            tag,
            _comm,
            &status );
}

template<typename T>
WxMsgStatus 
WxMpiMsgTmpl<T>::startRecv(unsigned num, unsigned sendRank, int tag)
{
  WxMpiMsgStatus_v *ms = new WxMpiMsgStatus_v();
  ms->data = new T[num];
  MPI_Irecv(ms->data,
            num,
            WxMpiTraits<T>::mpiType(),
            sendRank,
            tag,
            _comm,
            &ms->request);
  return ms;
}

template<typename T>
void
WxMpiMsgTmpl<T>::allReduce(unsigned num, T* sendBuff, T* recvBuff, WxMsgOp op)
{
  MPI_Allreduce(sendBuff,
                recvBuff,
                num,
                WxMpiTraits<T>::mpiType(),
                _ops.getOp(op),
                _comm);
}

// instantiate classes for all types define in wxdatatypes.h
template class WxMpiMsgTmpl<char>;
template class WxMpiMsgTmpl<unsigned char>;
template class WxMpiMsgTmpl<short>;
template class WxMpiMsgTmpl<unsigned short>;
template class WxMpiMsgTmpl<int>;
template class WxMpiMsgTmpl<unsigned int>;
template class WxMpiMsgTmpl<long>;
template class WxMpiMsgTmpl<unsigned long>;
template class WxMpiMsgTmpl<float>;
template class WxMpiMsgTmpl<double>;
template class WxMpiMsgTmpl<long double>;
template class WxMpiMsgTmpl<long long int>;
