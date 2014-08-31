#ifndef __wxlogstreamstrm__h__
#define __wxlogstreamstrm__h__

// WarpX lib includes

// std includes
#include <sstream>
#include <iostream>

// forward declare WxLogger class
class WxLoggerStream;
class WxLogger;

class WxLogStreamStrm {
  public:
    friend class WxLogStream;

/**
 * Output supplied value
 *
 * @param val value to output
 * @return reference to this stream object
 */
    template <typename T>
    WxLogStreamStrm& operator<<(T val) {
      _sstrm.str("");
      _sstrm << val;
      this->_logIt(_sstrm);
      return *this;
    }

/**
 * I/O for manipulators
 *
 * @param p manipulator object
 * @return reference to this stream object
 */
    WxLogStreamStrm&
    operator<<(std::ostream& (*p)(std::ostream&)) {
      _sstrm.str("");
      _sstrm << p;
      this->_logIt(_sstrm);
      return *this;
    }

/**
 * I/O for manipulators
 *
 * @param p manipulator object
 * @return reference to this stream object
 */
    WxLogStreamStrm& 
    operator<<(std::ios& (*p)(std::ios&)) {
      _sstrm.str("");
      _sstrm << p;
      this->_logIt(_sstrm);
      return *this;
    }

  private:
/** 
 * Ctor is private so only logg-stream can make instances.
 */
    WxLogStreamStrm(WxLogger* log, int level);

/**
 * Copy ctor is private to avoid copying
 */
    WxLogStreamStrm(const WxLogStreamStrm&);

    std::ostringstream _sstrm;
    WxLogger *_logger;
    int _level;

    void _logIt(const std::ostringstream& str);
};

#endif //  __wxlogstreamstrm__h__
