// WarpX lib includes
#include "wxobject.h"

WxObject::WxObject()
  : _name("noname"), _io(0), _msg(0)
{
}

WxObject::WxObject(const std::string& name)
  : _name(name), _io(0), _msg(0)
{
}

WxObject::~WxObject() 
{
}

void
WxObject::setIo(WxIoBase& io)
{
  _io = &io;
}

void
WxObject::setMsg(WxMsgBase& msg)
{
  _msg = &msg;
}

WxIoBase&
WxObject::getIo()
{
  return *_io;
}

WxMsgBase&
WxObject::getMsg()
{
  return *_msg;
}

void
WxObject::setup(const WxCryptSet& wxc) 
{
  _name = wxc.name();
}

void
WxObject::init(Vec X)
{ // by default do nothing
}

void
WxObject::load(WxIoBase& io, const WxIoNodeType& grpNode)
{ // by default do nothing
}

void
WxObject::dump(WxIoBase& io, WxIoNodeType& grpNode)
{ // by default do nothing
}

void
WxObject::finishBuild()
{ // by default do nothing
}

std::string
WxObject::name() const 
{
  return _name;
}

void 
WxObject::setName(const std::string& nm)
{
  _name = nm;
}
