#ifndef __wxlogger__h__
#define __wxlogger__h__

// WarpX includes
#include "wxlogrecordhandler.h"
#include "wxexcept.h"

// std includes
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

// forward declare WxLogStream class
class WxLogStream;

template<class T>
class WxLoggerBase
{
  public:
    typedef std::map<std::string, T*, std::less<std::string> > LoggerMap_t;
    typedef std::pair<std::string, T*> LoggerPair_t;

/**
 * Returns a logger with a given name. If such a logger does not exist
 * a new logger with that name is created and returned. Note that
 * calling this method with the same name will always return the same
 * logger object. Thus once loggers are created they become global
 * objects.
 */    
    static T* get(const std::string& nm) {

      if (!_loggers)
        _loggers = new LoggerMap_t();

      // check if logger with given name is in map
      typename LoggerMap_t::iterator i = _loggers->find(nm);
      if (i != _loggers->end())
        return i->second;

      // it is not, so start creation process
        
      // find location of '.'
      int len = 0;
      for (int i=nm.size(); i>=0; --i, ++len)
        if ( (nm[i]=='.') || (i==0) )
        {
          // create logger with given name
          T *l = new T(nm);
          _loggers->insert(LoggerPair_t(nm, l));
          // set its parent 
          if (i!=0)
            l->_parent = WxLoggerBase<T>::get(nm.substr(0, nm.size()-len));
          return l;
        }
      // this return statement is never executed
      return 0;
    }

/**
 * Returns a logger with a given name. If such a logger does not exist
 * an exception is thrown.
 */
    static T* getSafe(const std::string& nm) {

      if (!_loggers)
        _loggers = new LoggerMap_t();

      typename LoggerMap_t::iterator i = _loggers->find(nm);
      if (i != _loggers->end())
        return i->second;
      // thow exception
      WxExcept wxe;
      wxe << "Logger " << nm << " not found";
      throw wxe;
    }

/**
 * Delete all loggers registered in the system.
 */
    static void cleanUp() {
      typename LoggerMap_t::iterator itr;
      for (itr=_loggers->begin(); itr!=_loggers->end(); ++itr)
        delete itr->second;
    }

  private:
    static LoggerMap_t *_loggers;
};
// initialize maps
template<class T>
std::map<std::string, T*, std::less<std::string> >* WxLoggerBase<T>::_loggers = NULL;

class WxLogger : public WxLoggerBase<WxLogger>
{
  public:
    friend class WxLoggerBase<WxLogger>;

    enum eLevels {
      NOTSET = 0,
      DEBUG = 10,
      INFO = 20,
      WARNING = 30,
      ERROR = 40,
      CRITICAL = 50,
      DISABLED = 1000
    };

    // map type to map logger level to its string representation
    typedef std::map<unsigned, std::string, std::less<unsigned> > LevelMap_t;
    typedef std::pair<unsigned, std::string> LevelPair_t;

    typedef std::map<std::string, eLevels, std::less<std::string> > StringMap_t;
    typedef std::pair<std::string, eLevels>  StringPair_t;

/**
 * Creates a new logger with given name. The verbosity is set
 * NOTSET. However, by default no handlers are installed so no
 * messages will be logged.
 */
    WxLogger(const std::string &name);

/**
 * Delete logger and all log record handlers
 */
    virtual ~WxLogger();

/**
 * Log a debug message
 */
    void debug(const std::string& msg) const;

/**
 * Log a info message
 */
    void info(const std::string& msg) const;

/**
 * Log a warning message
 */
    void warning(const std::string& msg) const;

/**
 * Log a error message
 */
    void error(const std::string& msg) const;

/**
 * Log a critical message
 */
    void critical(const std::string& msg) const;

/**
 * Generic messsage logger.
 */
    void log(const std::string& msg, eLevels withLevel) const;

/**
 * Set verbosity level
 */
    void setLevel(enum eLevels level);

/**
 * Set verbosity level passing a string
 */
    void setLevel(const std::string& level);

/**
 * Get verbosity level
 */
    eLevels getLevel() const;

/**
 * Get verbosity level as a string
 */
    std::string getLevelStr();

/**
 * Add a new handler to the logger.
 *
 * @param handler Instance of class WxLogRecordHandler. Note that the
 * handler is not owned by the logger and hence is never deleted.
 */
    void addHandler(WxLogRecordHandler *handler);

/**
 * Disables all logging to this logger. Logging at the original level
 * can be resumed by calling the 'enable' method.
 */
    void disable();

/**
 * Re-enables logging to this logger. The verbosity level is set to
 * the one prior to the disable call.
 */
    void enable();

/**
 * Return stream to log debug messages
 */
    WxLogStream getDebugStream();

/**
 * Return stream to log info messages
 */
    WxLogStream getInfoStream();

/**
 * Return stream to log warning messages
 */
    WxLogStream getWarningStream();

/**
 * Return stream to log error messages
 */
    WxLogStream getErrorStream();

/**
 * Return stream to log critical messages
 */
    WxLogStream getCriticalStream();

  private:
    WxLogger *_parent; // parent logger
    std::string _name;
    eLevels _level, _oldLevel;
    LevelMap_t _levelMap;
    StringMap_t _stringMap;
    std::vector<WxLogRecordHandler*> _handlers;

};

#endif // __wxlogger__h__
