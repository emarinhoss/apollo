#ifndef __wxsingleton__
#define __wxsingleton__

template<class T>
struct WxCreateUsingNew
{
/**
 * Create object using the new operator
 *
 * @return pointer to newly creted object
 */
    static T* create() {
      return new T;
    }
};

/**
 * Provides a class to create singleton from supplied type
 */
template<class T,
         template <class> class CreationPolicy = WxCreateUsingNew
         >
class WxSingleton
{
  public:
/**
 * Returns reference to unique copy of the object of type T.
 *
 * @return Reference to object stored in singleton
 */
    static T& instance() {
      if (!_instance)
        _instance = CreationPolicy<T>::create();
      return *_instance;
    }

  private:
/** Global object stored as singleton */
    static T* _instance;
};
// initialize the instance to NULL
template<class T, template <class> class CreationPolicy> 
T* WxSingleton<T, CreationPolicy>::_instance = 0;

#endif //  __wxsingleton__
