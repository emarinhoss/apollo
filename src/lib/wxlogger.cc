// WarpX lib includes
#include "wxlogger.h"
#include "wxlogstream.h"

WxLogger::WxLogger(const std::string &name) 
  : _parent(0), _name(name), _level(NOTSET), _oldLevel(NOTSET) 
{

  // construct level -> string mapping
  _levelMap.insert( LevelPair_t(NOTSET, "notset") );
  _levelMap.insert( LevelPair_t(DEBUG, "debug") );
  _levelMap.insert( LevelPair_t(INFO, "info") );
  _levelMap.insert( LevelPair_t(WARNING, "warning") );
  _levelMap.insert( LevelPair_t(ERROR, "error") );
  _levelMap.insert( LevelPair_t(CRITICAL, "critical") );
  _levelMap.insert( LevelPair_t(DISABLED, "disabled") );

  // construct string -> level mapping
  _stringMap.insert( StringPair_t("notset", NOTSET) );
  _stringMap.insert( StringPair_t("debug", DEBUG) );
  _stringMap.insert( StringPair_t("info", INFO) );
  _stringMap.insert( StringPair_t("warning", WARNING) );
  _stringMap.insert( StringPair_t("error", ERROR) );
  _stringMap.insert( StringPair_t("critical", CRITICAL) );
  _stringMap.insert( StringPair_t("disabled", DISABLED) );
}

WxLogger::~WxLogger()
{
  std::vector<WxLogRecordHandler*>::iterator itr;
  for (itr=_handlers.begin(); itr!=_handlers.end(); ++itr)
    delete *itr;
}

void 
WxLogger::debug(const std::string& msg) const 
{
  log(msg, DEBUG);
}

void 
WxLogger::info(const std::string& msg) const 
{
  log(msg, INFO);
}

void 
WxLogger::warning(const std::string& msg) const 
{
  log(msg, WARNING);
}

void 
WxLogger::error(const std::string& msg) const 
{
  log(msg, ERROR);
}

void
WxLogger::critical(const std::string& msg) const 
{
  log(msg, CRITICAL);
}

void 
WxLogger::setLevel(enum eLevels level) 
{
  _oldLevel = _level = level;
}

void 
WxLogger::setLevel(const std::string& level) 
{
  StringMap_t::const_iterator i =
    _stringMap.find(level);
  // set level, defaulting to NOTSET if incorrect string passed
  _oldLevel = _level = 
    ((i != _stringMap.end()) ? (*i).second : NOTSET);
}

WxLogger::eLevels 
WxLogger::getLevel() const 
{
  return _level;
}

std::string 
WxLogger::getLevelStr() 
{
  return _levelMap[_level];
}

void 
WxLogger::addHandler(WxLogRecordHandler *handler) 
{
  _handlers.push_back(handler);
}

void 
WxLogger::disable() 
{
  _oldLevel = _level;
  _level = DISABLED;
}

void 
WxLogger::enable() 
{
  _level = _oldLevel;
}

WxLogStream 
WxLogger::getDebugStream() 
{
  return WxLogStream(this, DEBUG);
}

WxLogStream 
WxLogger::getInfoStream()
{
  return WxLogStream(this, INFO);
}

WxLogStream 
WxLogger::getWarningStream()
{
  return WxLogStream(this, WARNING);
}

WxLogStream 
WxLogger::getErrorStream()
{
  return WxLogStream(this, ERROR);
}

WxLogStream 
WxLogger::getCriticalStream()
{
  return WxLogStream(this, CRITICAL);
}

void
WxLogger::log(const std::string& msg, eLevels withLevel) const 
{
  // check if level of logger is sufficient to log this message
  if (_level <= withLevel)
  {
    // send message to each handler registered with logger
    std::vector<WxLogRecordHandler*>::const_iterator i;
    for (i=_handlers.begin(); i!=_handlers.end(); ++i) 
      (*i)->write(msg);
  }
  // now send this same message to the parent if there is one
  if (_parent)
    _parent->log(msg, withLevel);
}

