#ifndef __wxany__h__
#define __wxany__h__

// WarpX includes

// std includes
#include <algorithm>
#include <typeinfo>

/**
 * Class WxAny is based on the "any" class described in "Valued
 * Conversion", Kevlin Henney, C++ Report, July-August 2000, pages
 * 37--40. WxAny should be used for objects with simple types which
 * support a copy ctor.
 *
 * Example usage is:
 *
 * WxAny any(22); // store integer 22
 * cout << wx_any_cast<int>(any) << endl; // prints 22
 *
 *  any = WxAny(string("Hello World"); // store string
 *  cout << wx_any_cast<string>(any) << endl; // prints "Hello World"
 */
class WxAny
{
  public:

/**
 *  Create empty object
 */
    WxAny()
      : content(0) {
    }

/**
 * Create new WxAny from an object of a given type
 *
 *  @param value Value of object stored.
 */
    template<typename VALUETYPE>
    WxAny(const VALUETYPE& value)
      : content(new _WxHolder<VALUETYPE>(value)) {
    }

/**
 * Copy ctor
 */
    WxAny(const WxAny& other)
      : content(other.content ? other.content->clone() : 0) {
    }

    ~WxAny() {
      delete content;
    }

/**
 * Swap contents of object with contents of supplied object
 *
 * @param rhs Replace the value in this WxAny with 'rhs'
 */
    WxAny& swap(WxAny& rhs) {
      std::swap(content, rhs.content);
      return *this;
    }

/**
 * Assignment operator: use WxAny object to create a new object
 *
 * @param rhs WxAny to assign from
 */
    WxAny& operator=(const WxAny& rhs) {
      WxAny(rhs).swap(*this);
      return *this;
    }

/**
 * Assignment operator: use VALUETYPE object to create a new object
 */
    template<typename VALUETYPE>
    WxAny& operator=(const VALUETYPE& rhs) {
      WxAny(rhs).swap(*this);
      return *this;
    }

/**
 * Is object empty?
 */
    bool empty() const {
      return !content;
    }

/**
 * Type_info for this object's held value
 */
    const std::type_info& type() const {
      return content ? content->type() : typeid(void);
    }

/**
 * Convert the held object to a pointer. Do not use this function
 * directly without checking if the returned pointer is not NULL.
 */
    template<typename VALUETYPE>
    const VALUETYPE* to_ptr() const {
      // return pointer to held object or NULL if wrong type
      return 
        type() == typeid(VALUETYPE) 
        ? &static_cast<_WxHolder<VALUETYPE> *> (content)->held
        : 0
        ;
    }

/**
 * Convert help object to a void *. This is not a safe operation as
 * it breaks typechecking.
 */
    const void* to_void_ptr() const {
      return content->void_ptr();
    }

/**
 * Extract the data from the WxAny object. 
 *
 * If the type specified by the template 'VALUETYPE' is not the
 * correct type of the object stored a std::bad_cast exception is
 * thrown.
 */
    template<typename VALUETYPE>
    VALUETYPE to_value() const {
      // return object or throw exception if types do not match
      const VALUETYPE * result =
        type() == typeid(VALUETYPE) 
        ? &static_cast<_WxHolder<VALUETYPE> *> (content)->held
        : 0
        ;
      // if it is not null, return value else thow exception
      if (result)
        return *result;
      else
        throw std::bad_cast();  
    }

  private:

    class _WxPlaceHolder
    {
      public:
        // dtor
        virtual ~_WxPlaceHolder() {
        }

        // typeinfo: must be provided by children
        virtual const std::type_info& type() const = 0;

        // make a copy: must be provided by children
        virtual _WxPlaceHolder* clone() const = 0;

        // return object as a void *
        virtual const void * void_ptr() const = 0;
        
    };

    template<typename VALUETYPE>
    class _WxHolder : public _WxPlaceHolder
    {
      public:
        _WxHolder(const VALUETYPE& value)
          : held(value) {
        }

        // type_info about VALUETYPE
        virtual const std::type_info& type() const {
          return typeid(VALUETYPE);
        }

        // make a copy of held object
        virtual _WxPlaceHolder* clone() const {
          return new _WxHolder(held);
        }

        // return object as a void *
        virtual const void * void_ptr() const {
          return (const void*) &held;
        }

        VALUETYPE held; // held object 
    };

    _WxPlaceHolder *content;

};

/**
 * Extract the data from the WxAny object. 
 *
 * @param operand WxAny object from which to extract value. If the
 * type specified by the template 'VALUETYPE' is not the correct type
 * of the object stored a std::bad_cast exception is thrown.
 */
template<typename VALUETYPE>
VALUETYPE
wx_any_cast(const WxAny& operand)
{
  // return object or throw exception if types do not match
  // return object or throw exception if types do not match
  const VALUETYPE * result = operand.template to_ptr<VALUETYPE>();

  // if it is not null, return value else thow exception
  if (result)
    return *result;
  else
    throw std::bad_cast();
}

/**
 * Extract the data from the WxAny object. 
 *
 * @param operand WxAny object from which to extract value. If the
 * type specified by the template 'VALUETYPE' is not the correct type
 * of the object stored a std::bad_cast exception is thrown.
 */
template<typename VALUETYPE>
VALUETYPE
wx_any_cast(WxAny& operand)
{
  // return object or throw exception if types do not match
  const VALUETYPE * result = operand.template to_ptr<VALUETYPE>();

  // if it is not null, return value else thow exception
  if (result)
    return *result;
  else
    throw std::bad_cast();
}

#endif // __wxany__h__
