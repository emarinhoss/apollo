#ifndef __wxeventhandler__
#define __wxeventhandler__

// std includes
#include <string>
#include <vector>
#include <map>
#include <iostream>

// lib includes
#include <wxeventdispatchertmpl.h>

/**
 * Class to handle events of different types. The types of event
 * should be specified in the type list.
 */
template <typename CH, typename T>
class WxEventHandler : public WxEventHandlerTmpl<T>
{
  public:
/**
 * Destroy object
 */
    ~WxEventHandler()
    {
      // unregister ourself before dying
      disp->unregisterHandler(name, index);
    }

/**
 * Register callback function for an event.
 *
 * @param nm Name of event
 * @param ed Reference to event dispatcher object
 * @param obj Pointer to object containing callback
 * @param func Pointer to member function for callback
 */
    void registerCallBack(const std::string& nm, WxEventDispatcherTmpl<T>& ed, CH* obj, void(CH::*func)(const T&) ) 
    {
      name = nm;
      handlerObj = obj;
      handlerFunc = func;
      // register ourself with dispatcher
      index = ed.registerHandler(name, this);
      disp = &ed; // store pointer to dispatcher
    }

/**
 * Run callback method registered for this event
 *
 * @param event Event to dispatch to handler
 */
    void runCallBack(const T& event) 
    {
      ((handlerObj)->*(handlerFunc))(event);
    }

  private:
/** Index to registers handler */
    unsigned index;
/** Name of event */
    std::string name;
/** Pointer to handler object */
    CH* handlerObj;
/** Pointer to handler function */
    void (CH::*handlerFunc)(const T&);
/** Pointer to dispatcher object */
    WxEventDispatcherTmpl<T> *disp;
};

#endif //  __wxeventhandler__
