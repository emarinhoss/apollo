#ifndef __wxstepperstatus__
#define __wxstepperstatus__

// std includes
#include <limits>

template <typename REAL>
class WxStepperStatus
{
  public:
/**
 * Create status object.
 *
 * @param status true if status is success
 * @param suggest suggested time-step
 */
    WxStepperStatus(bool status=true, REAL suggest=0.0) {
      _status = status;
      _suggestedDt = (suggest == 0.0) ?
        std::numeric_limits<REAL>::max() : suggest;
    }

/**
 * Set status of stepper
 *
 * @param status Status of stepper
 */    
    void setStatus(bool status) {
      _status = status;
    }

/**
 * Get status of stepper
 *
 * @return status of stepper
 */    
    bool getStatus() const {
      return _status;
    }

/**
 * Set suggested time-step
 *
 * @param dt suggested time-step
 */
    void setSuggestedDt(REAL dt) {
      _suggestedDt = dt;
    }

/**
 * Get suggested time-step
 *
 * @return suggested time-step
 */
    REAL getSuggestedDt() const {
      return _suggestedDt;
    }

  private:
/** Status of stepper */
    bool _status;
/** Suggested time-step */
    REAL _suggestedDt;
};

#endif // __wxstepperstatus__
