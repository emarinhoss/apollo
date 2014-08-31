// WarpX lib includes
#include "wxlogstream.h"
#include "wxlogger.h"


WxLogStream::WxLogStream(WxLogger* log, int level)
  : useCount(new int(1))
{
  strm = new WxLogStreamStrm(log, level);
}

WxLogStream::WxLogStream(const WxLogStream& ls)
  : strm(ls.strm)
{
  ++*ls.useCount;
  useCount = ls.useCount;
}

WxLogStream::~WxLogStream()
{
  if (--*useCount == 0)
  {
    delete useCount;
    delete strm;
  }
}
