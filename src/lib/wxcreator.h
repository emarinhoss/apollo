/**
 * The set of classes in this file provide a mechanism for creating
 * self registering classes. If class D derived from class B then
 *
 * WxCreator<D,B> d("mydclass");
 *
 * registers a factory for objects of type D. After that doing
 *
 * B *p = WxCreatorMap<B>::getNew("mydclass");
 *
 * returns a new object of type D. This way different classes derived
 * from the same base B can be selected at run time based on the
 * derived class registered name.
 */
#ifndef __wxcreator__h__
#define __wxcreator__h__

// WarpX includes
#include "wxexcept.h"

// std includes
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

template<typename B> class WxCreatorBase;

template<class B>
class WxCreatorMapBase
{
  public:
    typedef std::map<std::string, WxCreatorBase<B>*, std::less<std::string> > CreatorMap_t;
    typedef std::pair<std::string, WxCreatorBase<B>*> CreatorPair_t;

/**
 * Add a new creator object into the list of available creators.
 *
 * @param nm Name of the creator
 * @param b Pointer to the creator base object
 */
    static void addCreator(const std::string& nm, WxCreatorBase<B>* b) {
      if (!_creators)
        _creators = new CreatorMap_t();
      _creators->insert(CreatorPair_t(nm, b));
    }

/**
 * Remove a creator from the list
 *
 * @param nm Name of creator to remove
 */
    static void removeCreator(const std::string& nm) {
      _creators->erase(nm);
    }

  protected:
    static CreatorMap_t *_creators;
};

template<class B>
class WxCreatorMap : public WxCreatorMapBase<B>
{
  public:
    typedef std::map<std::string, WxCreatorBase<B>*, std::less<std::string> > CreatorMap_t;

/**
 * Get a new object whose creator has the given name. The returned
 * object points to the base class.
 *
 *   @param nm Name of the creator.
 */
    static B* getNew(const std::string& nm) {
      if (WxCreatorMapBase<B>::_creators)
      {
        typename CreatorMap_t::iterator i = WxCreatorMapBase<B>::_creators->find(nm);
        if (i != WxCreatorMapBase<B>::_creators->end())
          return i->second->getNew();
      }
      WxExcept wxe;
      wxe << "Creator for class " << nm << " not found";
      throw wxe;
    }

/**
 * Get a list of registered names.
 */
    static std::vector<std::string> registeredNames() {
      std::vector<std::string> names;
      if (WxCreatorMapBase<B>::_creators)
      { // loop over each entry, extracting its name and shove it
        // into name list
        typename CreatorMap_t::const_iterator i;
        for (i=WxCreatorMapBase<B>::_creators->begin();
             i!=WxCreatorMapBase<B>::_creators->end(); ++i)
          names.push_back( (*i).first );
      }
      return names;
    }

/**
 * Checks if creator with given name is registered.
 */
    static bool has(const std::string& nm) {
      if (WxCreatorMapBase<B>::_creators)
      {
        typename CreatorMap_t::iterator i = WxCreatorMapBase<B>::_creators->find(nm);
        if (i != WxCreatorMapBase<B>::_creators->end())
          return true;
      }
      return false;
    }
};

template<class B>
class WxCreatorBase
{
  public:

    WxCreatorBase(const std::string& nm) 
      : _name(nm) {
      WxCreatorMapBase<B>::addCreator(_name, this);
    }

    virtual ~WxCreatorBase() {
      WxCreatorMapBase<B>::removeCreator(_name);
    }

    virtual B* getNew() = 0;

  private:
    std::string _name;
};

template<class D, class B>
class WxCreator : public WxCreatorBase<B>
{
  public:
    WxCreator(const std::string& name)
      : WxCreatorBase<B>(name) {
    }

    B* getNew() {
      D *d = new D;
      return dynamic_cast<B*>(d);
    }
};

// initialize map
template <class B>
std::map<std::string, WxCreatorBase<B>*, std::less<std::string> >* WxCreatorMapBase<B>
::_creators = NULL;

#endif // __wxcreator__h__
