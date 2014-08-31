// WarpX lib includes
#include "wxlogstreamstrm.h"
#include "wxlogger.h"

WxLogStreamStrm::WxLogStreamStrm(WxLogger* log, int level)
  : _logger(log), _level(level) 
{
}

void
WxLogStreamStrm::_logIt(const std::ostringstream& str)
{
  _logger->log(str.str(), (WxLogger::eLevels)_level);
}

