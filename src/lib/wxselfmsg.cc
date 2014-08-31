// lib includes
#include "wxselfmsg.h"

WxSelfMsg::WxSelfMsg()
{
    // add templated messengers to the base class
    this->addMsg<char>
        ( new WxSelfMsgTmpl<char>() );
    this->addMsg<unsigned char>
        ( new WxSelfMsgTmpl<unsigned char>() );
    this->addMsg<short>
        ( new WxSelfMsgTmpl<short>() );
    this->addMsg<unsigned short>
        ( new WxSelfMsgTmpl<unsigned short>() );
    this->addMsg<int>
        ( new WxSelfMsgTmpl<int>() );
    this->addMsg<unsigned int>
        ( new WxSelfMsgTmpl<unsigned int>() );
    this->addMsg<long>
        ( new WxSelfMsgTmpl<long>() );
    this->addMsg<unsigned long>
        ( new WxSelfMsgTmpl<unsigned long>() );
    this->addMsg<float>
        ( new WxSelfMsgTmpl<float>() );
    this->addMsg<double>
        ( new WxSelfMsgTmpl<double>() );
    this->addMsg<long double>
        ( new WxSelfMsgTmpl<long double>() );
    this->addMsg<long long int>
        ( new WxSelfMsgTmpl<long long int>() );
}

WxMsgBase*
WxSelfMsg::createSubComm(const std::vector<int>& ranks) 
{
    return new WxSelfMsg();
}

void * 
WxSelfMsg::finishRecv(WxMsgStatus wxms)
{
    return 0;
}

bool
WxSelfMsg::checkRecv(WxMsgStatus wxms)
{
    return true;
}
