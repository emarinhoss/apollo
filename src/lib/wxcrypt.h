#ifndef __wxcrypt__h__
#define __wxcrypt__h__

// WarpX include
#include "wxany.h"
#include "wxdatatypes.h"
#include "wxexcept.h"
#include "wxtypelist.h"

// std includes
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

/**
 * WxCrypt provides a container to store/retrive name-value pairs. The
 * names are of type std::string while the associated values can be of
 * arbitrary types supporting the copy constructor. Note that for
 * efficiency reasons WxCrypt should be used to store only
 * "light-weight" objects like native-type objects (int, double,
 * string etc.) or small objects which are not expensive to copy.
 *
 * WxCrypt objects are immutable: values once added can not be removed
 * or modifed.
 */
class WxCrypt
{
  public:
    // types of string->WxAny map and pairs
    typedef std::map<std::string, WxAny, std::less<std::string> > AnyMap_t;
    typedef std::pair<std::string, WxAny> AnyPair_t;

    WxCrypt();
    virtual ~WxCrypt();

/**
 * Copy ctor: make a deep copy of 'crypt'
 */
    WxCrypt(const WxCrypt& crypt);

/**
 * Assignment operator: LHS is a deep copy of 'rhs'
 */
    WxCrypt& operator=(const WxCrypt& rhs);

/**
 * Insert a name, value pair. Returns true if the insertion worked,
 * false otherwise
 *
 * @param name Name of object
 * @param value Value of object
 */
    template<typename VALUETYPE>
    bool add(const std::string& name, VALUETYPE value) {
      addName<VALUETYPE>(name);
      return _values.insert( AnyPair_t(name, value) ).second;
    }

/**
 * Returns true if the 'name' exists in the crypt.
 *
 * @param name Name of object to check.
 */
    bool has(const std::string& name) const;

/**
 * Retrieve the value associated with name
 *
 * @param name Return value associated with 'name'
 */
    template <typename VALUETYPE>
    VALUETYPE get(const std::string& name) const {
      AnyMap_t::const_iterator i = _values.find(name);
      if (i != _values.end())
        return wx_any_cast<VALUETYPE>((*i).second);
      // value not found: throw an exception
      WxExcept wxe;
      wxe << "Value " << name << " not found";
      throw wxe;
    }

/**
 * Retrieve list of values associated with name
 *
 * @param name Return list of values associated with 'name'
 */
    template <typename VALUETYPE>
    std::vector<VALUETYPE> getVec(const std::string& name) const {
      // first fetch vector of WxAnys
      std::vector<WxAny> vals = 
        this->template get<std::vector<WxAny> >(name);
      // now convert it into vector of VALUETYPE
      std::vector<VALUETYPE> res;
      for (unsigned i=0; i<vals.size(); ++i)
        res.push_back( wx_any_cast<VALUETYPE>(vals[i]) );
      return res;
    }

/**
 * Get names of inserted types 
 *
 * @return list of names
 */
    template <typename VALUETYPE>
    std::vector<std::string> getNames() const {
      return wxTypeMapExtract<VALUETYPE, WxTypeToNames>(_typeToNames).names;
    }

  private:

    // container which maps types -> names
    template <typename T>
    struct WxTypeContainer {
        std::vector<std::string> names;
    };
    typedef WxTypeMap<WxDataTypes_t, WxTypeContainer> WxTypeToNames;

    WxTypeToNames _typeToNames;
    AnyMap_t _values;

/**
 * Add a name to type->names map for supplied type
 */
    template <typename VALUETYPE>
    void addName(const std::string& name) {
      wxTypeMapExtract<VALUETYPE>(_typeToNames).names.push_back(name);
    }

/**
 * Copy names from crypt type map to this object
 */
    template <typename T>
    void copyNames(const WxCrypt& crypt) {
      std::vector<std::string> nms = crypt.getNames<T>();
      std::vector<std::string>::const_iterator i;
      for (i=nms.begin(); i!=nms.end(); ++i)
        wxTypeMapExtract<T>(_typeToNames).names.push_back(*i);
    }
};

#endif // __wxcrypt__h__
