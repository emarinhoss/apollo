#ifndef __wxtimer__
#define __wxtimer__

// std includes
#include <iostream>
#include <string>
#include <sstream>
#include <ctime>

#ifdef _DO_USE_MPI_
#include <mpi.h>
#endif

/**
 * Timer class which keeps accurate time and also works for more than
 * 24 hour time frame.
 */
class WxTimer {
  public:
/**
 * Constructor
 */
    WxTimer()
      : start(0), end(0) {
    }

/**
 * Start the timer.
 */
    void startTimer() {
      startCoarse = time(0);
      start = clock();
    }

/**
 * Stop the timer.
 */
    void stopTimer() {
      end = clock();
      endCoarse = time(0);
    }

/**
 * Return number of days elapsed
 *
 * @return Days elapsed
 */
    unsigned daysElapsed() const {
      return (int) difftime(endCoarse, startCoarse)/(24*60*60);
    }

/**
 * Return number of seconds elapsed. This is always less than a day
 * and needs to be added to the number of days elasped.
 *
 * @return Number of seconds
 */
    double secondsElapsed() const {
      unsigned days = daysElapsed();
      if (days > 0) {
        double t1 = 24*60*60 - (double)start/CLOCKS_PER_SEC;
        double t2 = (double)end/CLOCKS_PER_SEC;
        return t1+t2;
      }
      else {
        return ((double)(end-start))/CLOCKS_PER_SEC;
      }
    }

/**
 * Return time elapsed as a string.
 *
 * @return String representing time elapsed
 */
    std::string timeElapsedAsString() const {
      std::ostringstream te;
      unsigned days = daysElapsed();
      double seconds = secondsElapsed();
      if (days >01) {
        if (days == 1)
          te << " 1 day and " << seconds << " seconds";
        else
          te << days << " days and " << seconds << " seconds";
      }
      else {
        te << seconds << " seconds";
      }
      return te.str();
    }

  private:
/** Start and end times. These wrap around each 24 hours */
    clock_t start, end;
/** Start and end times. These are only accurate to the closest second */
    time_t startCoarse, endCoarse;
};

#endif // __wxtimer__

