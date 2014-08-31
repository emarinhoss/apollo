#ifndef __wxeventdispathcertmpl__
#define __wxeventdispathcertmpl__

// lib includes
#include <wxtypelist.h>
#include <wxexcept.h>
#include <wxeventhandlertmpl.h>

// std includes
#include <map>
#include <vector>
#include <string>

template <typename T>
class WxEventDispatcherTmpl
{
  public:
/**
 * Register an event of given type T and name. Derived classes should
 * call this method to register each event they will dispatch. All
 * registration must be done before any distaches are made.
 *
 * @param T type of event
 * @param name Name of event
 */
    void registerEvent(const std::string& name) 
    {
      // fetch reference to list of handlers for event
      typename std::map<std::string, HandlerList>::iterator 
        itr = handlers.find(name);
      if (itr == handlers.end())
        // no such event: set empty list of handlers
        handlers[name] = HandlerList();
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
    unsigned registerHandler(const std::string& name, WxEventHandlerTmpl<T>* handler) 
    {
      // fetch reference to list of handlers for event
      typename std::map<std::string, HandlerList>::iterator 
        itr = handlers.find(name);
      if (itr == handlers.end())
      { // no such event
        WxExcept wxe(L"WxEventDispatcherTmpl::registerHandler : No such event");
        throw wxe;
      }
      // determine where this handler will be put
      unsigned idx = itr->second.handlerObjList.size();
      // add handler
      itr->second.handlerObjList.push_back(handler);
      return idx;
    }

/**
 * Unregister a handler at a given index location for an event
 *
 * @param name Name of event
 * @param idx Index of handler to unregister
 */
    void unregisterHandler(const std::string& name, unsigned idx)
    {
      // fetch reference to list of handlers for event
      typename std::map<std::string, HandlerList>::iterator 
        itr = handlers.find(name);
      // delete specified handler
      itr->second.handlerObjList.erase(
          itr->second.handlerObjList.begin() + idx);
    }

/**
 * Dispatch a event of given name.
 *
 * @param name Name of event to dispatch
 * @param event Event object to send to handlers
 */
    void dispatchEvent(const std::string& name, const T& event) 
    {
      // fetch reference to list of handlers for event
      typename std::map<std::string, HandlerList>::iterator 
        itr = handlers.find(name);
      if (itr == handlers.end())
      { // no such event
        WxExcept wxe(L"WxEventDispatcherTmpl::dispatchEvent : No such event");
        throw wxe;
      }
      // run callback methods for each handler
      for (unsigned i=0; i<itr->second.handlerObjList.size(); ++i)
        itr->second.handlerObjList[i]->runCallBack(event);
    }

  private:
/**
 * Stores a list of object/function-pointer pair for an event of given
 * type
 */
    struct HandlerList
    {
/** List of handler objects */
        std::vector<WxEventHandlerTmpl<T>* > handlerObjList;
    };

/** List of handlers associated with events */
    std::map<std::string, HandlerList > handlers;
};

#endif // __wxeventdispathcertmpl__
