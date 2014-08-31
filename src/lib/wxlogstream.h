#ifndef __wxlogstream__h__
#define __wxlogstream__h__

// WarpX lib includes
#include <wxlogstreamstrm.h>

// std includes
#include <sstream>
#include <iostream>

class WxLogStream {
  public:
/** Create new stream object **/
    WxLogStream(WxLogger* log, int level);

/** Copy ctor */
    WxLogStream(const WxLogStream& ls);

/** Delete stream */
    ~WxLogStream();

/**
 * Output supplied value
 *
 * @param val value to output
 * @return reference to this stream object
 */
    template <typename T>
    WxLogStream& operator<<(T val) {
      strm->operator<<(val);
      return *this;
    }

/**
 * I/O for manipulators
 *
 * @param p manipulator object
 * @return reference to this stream object
 */
    WxLogStream&
    operator<<(std::ostream& (*p)(std::ostream&)) {
      strm->operator<<(p);
      return *this;
    }

/**
 * I/O for manipulators
 *
 * @param p manipulator object
 * @return reference to this stream object
 */
    WxLogStream& 
    operator<<(std::ios& (*p)(std::ios&)) {
      strm->operator<<(p);
      return *this;
    }

  private:
/** For reference counting **/
    mutable int *useCount;
/** Pointer to stream which does output */
    WxLogStreamStrm *strm;
};

#endif //  __wxlogstream__h__
