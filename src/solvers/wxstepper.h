#ifndef __wxstepper__
#define __wxstepper__

// WarpX lib includes
#include <wxobject.h>
#include <petsc.h>

// WarpX solver includes
#include <wxstepperstatus.h>

/**
 * Base class for objects which can be advanced in time. Classes
 * derived from this must provide the step(dt) method which advances
 * the state of the object by the given time step dt. Objects can
 * obtain the current time at which they are called by using the
 * getCurrentTime() method.
 */
template <typename REAL>
class WxStepper : public WxObject
{
  public:
/** 
 * Create a new stepper with given name 
 *
 * @param name Name of stepper
 */
    WxStepper(const std::string& name);

/** Dtor: destroy object */
    virtual ~WxStepper();

/**
 * Set time-step for the step method.
 *
 * @param dt Time-step
 */
    void setDt(REAL dt);

/**
 * Get time-step for the step method.
 *
 * @return Time-step
 */
    REAL getDt() const {
      return _dt;
    }

/**
 * Set current time.
 *
 * @param tcurr Current time
 */
    void setCurrentTime(REAL tcurr);

/**
 * Get the current time
 *
 * @return Current time
 */
    REAL getCurrentTime() const {
      return _currTime;
    }

/**
 * Advance the stepper by given time step. If this step failed and if
 * the parent solver is running in variable time-stepping mode, it
 * will be called again with a smaller time-step. A time-step
 * suggestion must be provided. The status of the step should be
 * encoded in the return object. See WxStepperStatus for details on
 * what fields can be set.
 *
 * @param dt time step to advance by
 * @return Status of stepper
 */
    virtual WxStepperStatus<REAL> step(REAL dt, Vec in, Vec out)=0;

  private:
/** Current time at which stepper is called */
    REAL _currTime;
/** Current time-step */
    REAL _dt;
};

#endif // __wxstepper__
