#ifndef __wxcryptset__h__
#define __wxcryptset__h__

// WarpX includes
#include "wxcrypt.h"
#include "wxcryptsetlexer.h"

// std incudes
#include <istream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>
#include <cstring>

/**
 * WxCryptSet extends WxCrypt by providing, in addition to name-value
 * pairs, an set of named WxCryptSets, thus providing a powerful way
 * of representing hierarchical data. In effect the WxCryptSet is a
 * n-ary tree with non-unique keys.
 *
 * Copying (or assignment of) a WxCryptSet can be expensive as the
 * complete tree is copied (deep-copy semantics). In general this is
 * not a problem as the creation of the set is usually done at
 * start-up time and the resulting set is never modified.
 *
 * Like its parent WxCrypt, WxCryptSet is also immutable: sets once
 * added can not be deleted or modified.
 */
class WxCryptSet : public WxCrypt
{
  public:
    // types of string->WxCryptSet map and pairs
    typedef std::map<std::string, WxCryptSet*, std::less<std::string> > CryptSetMap_t;
    typedef std::pair<std::string, WxCryptSet*> CryptSetPair_t;

    // types of type -> list of WxCryptSets name
    typedef std::map<std::string, std::vector<std::string> > TypeMap_t;
    typedef std::pair<std::string, std::vector<std::string> > TypePair_t;

/**
 * Create nameless cryptset
 */
    WxCryptSet();

/**
 * Create empty cryptset with given name.
 *
 * @param name Name of the set
 */
    WxCryptSet(const std::string& name);

/**
 * Create cryptset by parsing input from specified stream
 *
 * @param istr The input stream from which the cryptset is constructed
 */
    WxCryptSet(std::istream& istr);

    virtual ~WxCryptSet();

/**
 * Copy ctor: makes a deep copy of the supplied set
 */
    WxCryptSet(const WxCryptSet& wxc);

/**
 * Assignment operator: makes a deep copy of supplied set
 */
    WxCryptSet& operator=(const WxCryptSet& wxc);

/**
 * Name of crypt set
 */
    std::string name() const {
      return _name;
    }

/**
 * Add a new WxCryptSet.
 *
 * @param wxc Add given cryptset to this set.
 */
    void addSet(const WxCryptSet& wxc);

/**
 * Returns true if specified cyrptset exists
 *
 * @param name Name of cryptset to test for.
 */
    bool hasSet(const std::string& name) const;

/**
 * Get cryptset with given name. The returned set is immutable.
 *
 * @param name Name of cryptset to get.
 */
    const WxCryptSet& getSet(const std::string& name) const;

/**
 * Return list of cryptset names with the given type.
 *
 * @param type Type name
 * @return list of cryptset of the supplied type
 */
    std::vector<std::string> getNamesOfType(const std::string& type) const;

  private:
    std::string _name;
    CryptSetMap_t _cryptSets;
    TypeMap_t _typeMap;

/**
 * Change name of set
 *
 * @param name Name of cryptset
 */
    void _setName(const std::string& name) {
      _name = name;
    }

    // class to parse input stream and initialize cryptset.
    class _WxParser
    {
      public:
        _WxParser(std::istream& istr, WxCryptSet *wxcs);
        _WxParser(_WxParser *wxl, WxCryptSet *wxcs);

        // symbol lookup functions
        bool accept(int sym);
        bool expect(int sym);

        // non-terminals 
        void input_file();
        void wx_decl();
        void wx_decl_top();
        void wx_decl_end();
        void wx_decl_list();
        void wx_decl_list_rest();
        void wx_decl_value();
        void wx_value();
        void wx_value_list();
        void wx_value_sub_list();
        void wx_value_sub_list_rest();

        WxCryptSetLexer _fl; // lexer
        WxCryptSet* _wxcs; // pointer to cryptset being modified
        int _sym; // current symbol read
        std::string _tokenStr; // char pointer to token string
        std::vector<WxAny> _tokens; // list of tokens
        std::vector<WxAny> _valueList; // parse type of list
        std::string _idName; // name of attribute being parsed
        bool isParsingList; // are we parsing a list?
    };

/**
 * Parse input stream into cryptset
 */
    void _parse(_WxParser& fl);
};

#endif // __wxcryptset__h__
