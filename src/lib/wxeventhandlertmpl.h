#ifndef __wxeventhandlertmpl__
#define __wxeventhandlertmpl__

/**
 * Provides a place holder for classes that want to handle an event of
 * type T.
 */
template <typename T>
class WxEventHandlerTmpl
{
  public:
/** Destroy object */
    virtual ~WxEventHandlerTmpl() 
    {
    }

    virtual void runCallBack(const T& event) 
    { // defer to WxEventHandler class      
    }
};

#endif // __wxeventhandlertmpl__
