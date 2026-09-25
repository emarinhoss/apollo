// WarpX lib includes
#include "wxtypelist.h"
#include "wxdatatypes.h"
#include "wxselfmsgtmpl.h"

template<typename T>
void 
WxSelfMsgTmpl<T>::send(const std::vector<T>& arr, unsigned recvRank, int tag) 
{
}

template<typename T>
void 
WxSelfMsgTmpl<T>::send(unsigned num, T* arr, unsigned recvRank, int tag) 
{
}

template<typename T>
void 
WxSelfMsgTmpl<T>::recv(unsigned num, std::vector<T>& array, unsigned sendRank, int tag) 
{
}

template<typename T>
void 
WxSelfMsgTmpl<T>::recv(unsigned num, T* array, unsigned sendRank, int tag) 
{
}

template<typename T>
WxMsgStatus 
WxSelfMsgTmpl<T>::startRecv(unsigned num, unsigned sendRank, int tag)
{
  return new WxSelfMsgStatus_v();
}

template<typename T>
void
WxSelfMsgTmpl<T>::allReduce(unsigned num, T* sendBuff, T* recvBuff, WxMsgOp op)
{
  for (unsigned i=0; i<num; ++i)
    recvBuff[i] = sendBuff[i];
}

// instantiate classes for all types define in wxdatatypes.h
template class WxSelfMsgTmpl<char>;
template class WxSelfMsgTmpl<unsigned char>;
template class WxSelfMsgTmpl<short>;
template class WxSelfMsgTmpl<unsigned short>;
template class WxSelfMsgTmpl<int>;
template class WxSelfMsgTmpl<unsigned int>;
template class WxSelfMsgTmpl<long>;
template class WxSelfMsgTmpl<unsigned long>;
//template class WxSelfMsgTmpl<float>;
template class WxSelfMsgTmpl<double>;
template class WxSelfMsgTmpl<long double>;
template class WxSelfMsgTmpl<long long int>;
