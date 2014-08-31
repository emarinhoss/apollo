#ifndef __wxeventdispathcer__
#define __wxeventdispathcer__

// lib includes
#include <wxtypelist.h>
#include <wxeventdispatchertmpl.h>

// std includes
#include <map>
#include <vector>
#include <string>

/**
 * Base class for classes which wish to dispatch events of different
 * types. The types of event should be specified in the type list.
 */
template <typename TYPELIST>
class WxEventDispatcher : public WxTypeMap<TYPELIST, WxEventDispatcherTmpl>
{
  public:    
/**
 * Create a new dispatcher object
 */
    WxEventDispatcher() 
    {
      // this is delibrately empty
    }

/**
 * Register an event of given type T and name. Derived classes should
 * call this method to register each event they will dispatch. All
 * registration must be done before any distaches are made.
 *
 * @param T type of event
 * @param name Name of event
 */
    template <typename T>
    void registerEvent(const std::string& name) 
    { // redirect to parent
      WxEventDispatcherTmpl<T>::registerEvent(name);
    }

/**
 * Register a handler for the an event dispatched by this
 * dispatcher. This method allows objects to register functions which
 * will be called when the event they want to handle are dispatched.
 *
 * @param name Name of event
 * @param handler Pointer to handler object
 * @return integer for use in un-registering the handler
 */
    template <typename T>
    unsigned registerHandler(const std::string& name, WxEventHandlerTmpl<T>* handler) 
    { // redirect to parent
      return WxEventDispatcherTmpl<T>::registerHandler(name, handler);
    }

/**
 * Dispatch a event of given name.
 *
 * @param name Name of event to dispatch
 * @param event Event object to send to handlers
 */
    template <typename T>
    void dispatchEvent(const std::string& name, const T& event) 
    { // redirect to parent
      WxEventDispatcherTmpl<T>::dispatchEvent(name, event);
    }

  private:
/**
 * Private copy ctor to prevent copying
 *
 * @param ed Event dispatcher
 */
    WxEventDispatcher(const WxEventDispatcher<TYPELIST>& ed);

/**
 * Private assignment operator to prevent assignment
 *
 * @param ed Event dispatcher
 */
    WxEventDispatcher& operator=(const WxEventDispatcher<TYPELIST>& ed);
};

#endif // __wxeventdispathcer__
