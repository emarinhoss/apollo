// WarpX includes
#include "wxcryptset.h"
#include "wxlogger.h"
#include "wxlogrecordhandler.h"
#include "wxexcept.h"

#include <sstream>

WxCryptSet::WxCryptSet()
  : WxCrypt(), _name("")
{}

WxCryptSet::WxCryptSet(const std::string& name)
  : WxCrypt(), _name(name)
{}

WxCryptSet::WxCryptSet(std::istream& istr)
  : WxCrypt(), _name("")
{
  _WxParser wxl(istr, this);
  _parse(wxl);
}

WxCryptSet::~WxCryptSet()
{
  CryptSetMap_t::iterator i;
  for (i=_cryptSets.begin(); i!=_cryptSets.end(); ++i)
    delete i->second;
  _cryptSets.erase( _cryptSets.begin(), _cryptSets.end() );
}

WxCryptSet::WxCryptSet(const WxCryptSet& wxc)
  : WxCrypt(wxc)
{
  CryptSetMap_t::const_iterator i;
  _name = wxc._name;

  for (i=wxc._cryptSets.begin(); i!=wxc._cryptSets.end(); ++i)
  {
    WxCryptSet *wxcp = new WxCryptSet(*(i->second));
    _cryptSets.insert( CryptSetPair_t( i->first, wxcp ) );
    // check if cryptset has Type attribute
    if (wxcp->has("Type"))
    {
      std::string type = wxcp->get<std::string>("Type");
      // add name to type map
      TypeMap_t::iterator i;
      i = _typeMap.find(type);
      if (i != _typeMap.end())
        (i->second).push_back(wxcp->name());
      else
      {
        std::vector<std::string> name;
        name.push_back(wxcp->name());
        _typeMap.insert( TypePair_t(type, name) );
      }
    }
  }
}

WxCryptSet&
WxCryptSet::operator=(const WxCryptSet& wxc)
{
  if (this==&wxc) return *this;

  // call base class assignment operator
  WxCrypt::operator= (wxc);

  CryptSetMap_t::const_iterator i;
  // delete all entries in map
  for (i=_cryptSets.begin(); i!=_cryptSets.end(); ++i)
    delete i->second;
  _cryptSets.erase( _cryptSets.begin(), _cryptSets.end() );

  // copy entries from wxc
  _name = wxc._name;
  for (i=wxc._cryptSets.begin(); i!=wxc._cryptSets.end(); ++i)
  {
    WxCryptSet *wxcp = new WxCryptSet(*(i->second));
    _cryptSets.insert( CryptSetPair_t( i->first, wxcp ) );
    // check if cryptset has Type attribute
    if (wxcp->has("Type"))
    {
      std::string type = wxcp->get<std::string>("Type");
      // add name to type map
      TypeMap_t::iterator i;
      i = _typeMap.find(type);
      if (i != _typeMap.end())
        (i->second).push_back(wxcp->name());
      else
      {
        std::vector<std::string> name;
        name.push_back(wxcp->name());
        _typeMap.insert( TypePair_t(type, name) );
      }
    }
  }

  return *this;
}

void 
WxCryptSet::addSet(const WxCryptSet& wxc) 
{
  // make a deep copy of the set
  WxCryptSet *wxcp = new WxCryptSet(wxc);
  _cryptSets.insert( CryptSetPair_t(wxc.name(), wxcp) );
  // check if cryptset has Type attribute
  if (wxc.has("Type"))
  {
    std::string type = wxc.get<std::string>("Type");
    // add name to type map
    TypeMap_t::iterator i;
    i = _typeMap.find(type);
    if (i != _typeMap.end())
      (i->second).push_back(wxc.name());
    else
    {
      std::vector<std::string> name;
      name.push_back(wxc.name());
      _typeMap.insert( TypePair_t(type, name) );
    }
  }
}

bool
WxCryptSet::hasSet(const std::string& name) const 
{
  // retrieve the set with given name
  CryptSetMap_t::const_iterator i;
  i = _cryptSets.find(name);
  return (i != _cryptSets.end()) ? true : false;
}

const 
WxCryptSet& 
WxCryptSet::getSet(const std::string& name) const 
{
  // retrieve the set with given name
  CryptSetMap_t::const_iterator i;
  i = _cryptSets.find(name);
  if (i != _cryptSets.end())
    return *(i->second);
  WxExcept wxe;
  wxe << "Crypt set " << name << " not found";
  throw wxe;
}

std::vector<std::string> 
WxCryptSet::getNamesOfType(const std::string& type) const
{
  TypeMap_t::const_iterator i;
  i = _typeMap.find(type);
  if (i != _typeMap.end())
    return i->second;
  else
    return std::vector<std::string>();
}

//
// Parser code: a simple recursive decent parser implementing LL(1)
// grammar in file inpgrammar_ll.
//

WxCryptSet::_WxParser::_WxParser(std::istream& istr, WxCryptSet *wxcs)
  : _fl(istr),  _wxcs(wxcs)
{
  // initialize parsing by reading first token
  _sym = _fl.yylex();
  _tokenStr = _fl.YYText();
}

WxCryptSet::_WxParser::_WxParser(_WxParser *wxl, WxCryptSet *wxcs)
  : _fl(wxl->_fl), _wxcs(wxcs), _sym(wxl->_sym),
    _tokenStr(wxl->_tokenStr) 
{
}

bool
WxCryptSet::_WxParser::accept(int sym)
{
  if (sym == _sym)
  {
    // copy text of token 'sym' and advance the token stream
    _tokenStr = _fl.YYText();
    _sym = _fl.yylex();
    return true;
  }
  return false;
}

bool
WxCryptSet::_WxParser::expect(int sym)
{
  if ( accept(sym) )
    return true;

  std::ostringstream ss;
  ss << "Error: Unexpected symbol " << _tokenStr
     << "found  on line no "
     << _fl.lineno();
  WxLogger::get("wxcryptset-logger")->error(ss.str());
  throw WxExcept(ss.str());
}

void
WxCryptSet::_WxParser::wx_decl_top()
{
  expect(WX_LEFT_ANGLE);
  expect(WX_ID);
  _tokens.push_back( WxAny(_tokenStr) );
  WxLogger::get("wxcryptset-logger")->info(
      std::string("<") +
      _tokenStr +
      std::string(">\n"));
  expect(WX_RIGHT_ANGLE);
}

void
WxCryptSet::_WxParser::wx_decl_end()
{
  expect(WX_LEFT_ANGLE_FRONT_SLASH);
  expect(WX_ID);
  WxLogger::get("wxcryptset-logger")->info(
      std::string("</") +
      _tokenStr +
      std::string(">\n"));
  expect(WX_RIGHT_ANGLE);
}

void
WxCryptSet::_WxParser::wx_value()
{
  WxLogger *l = WxLogger::get("wxcryptset-logger");
  std::ostringstream ss;
  if (accept(WX_INT))
  {
    _tokens.push_back( WxAny(_fl.integer()) );
    if (!isParsingList)
      _wxcs->add<int>(_idName, _fl.integer());
    ss << "<int>" << _fl.integer() << "</int>\n";
    l->info(ss.str());
  }
  else if (accept(WX_REAL))
  {
    _tokens.push_back( WxAny(_fl.real()) );
    if (!isParsingList)
      _wxcs->add<double>(_idName, _fl.real());
    ss << "<real>" << _fl.real() << "</real>\n";
    l->info(ss.str());
  }
  else if (accept(WX_ID))
  {
    _tokens.push_back( WxAny(_tokenStr) );
    if (!isParsingList)
      _wxcs->add<std::string>(_idName, _tokenStr);
    l->info(std::string("<id>") + 
      _tokenStr + 
      std::string("</id>\n"));
  }
  else if (accept(WX_VALUE))
  {
    _tokens.push_back( WxAny(_tokenStr) );
    if (!isParsingList)
      _wxcs->add<std::string>(_idName, _tokenStr);
    l->info(std::string("<value>") + 
      _tokenStr + 
      std::string("</value>\n"));
  }
  else if (accept(WX_STRING))
  {
    _tokens.push_back( WxAny(_tokenStr) );
    if (!isParsingList)
      _wxcs->add<std::string>(_idName, _tokenStr);
    l->info(std::string("<string>") + 
      _tokenStr + 
      std::string("</string>\n"));
  }
  else
  {
    std::ostringstream ss;
    ss << "Error: Unknown token " << _tokenStr << " encountered on line no "
       << _fl.lineno() << std::endl;
    l->error(ss.str());
    WxExcept wxe;
    wxe << ss.str();
    throw wxe;
  }
}

void
WxCryptSet::_WxParser::wx_value_sub_list_rest()
{
  if ( accept(WX_COMMA) )
  {
    wx_value();
    // add value ...
    _valueList.push_back( _tokens.back() );
    _tokens.pop_back(); // .. and remove it
    wx_value_sub_list_rest();
  }
}

void
WxCryptSet::_WxParser::wx_value_sub_list()
{
  wx_value();
  // add value ...
  _valueList.push_back( _tokens.back() );
  _tokens.pop_back(); // .. and remove it
  wx_value_sub_list_rest();
}

void
WxCryptSet::_WxParser::wx_value_list()
{
  if (accept(WX_LEFT_BOX))
  {
    isParsingList = true;
    WxLogger::get("wxcryptset-logger")->info("<list>");
    // clear _valueList to insert list
    _valueList.erase( _valueList.begin(), _valueList.end() );
    wx_value_sub_list();
    expect(WX_RIGHT_BOX);
    WxLogger::get("wxcryptset-logger")->info("</list>\n");
    // add it to parsed data
    _tokens.push_back( _valueList );
    _wxcs->add(_idName, _valueList);
    isParsingList = false;
  }
  else
  {
    isParsingList = false;
    wx_value();
  }
}

void
WxCryptSet::_WxParser::wx_decl_value()
{
  if (accept(WX_ID))
  {
    std::string name = _tokenStr;
    _idName = name;
    expect(WX_EQUAL);
    WxLogger::get("wxcryptset-logger")->info(
        std::string("<") + name + std::string(">\n"));
    // get value
    wx_value_list();
    WxLogger::get("wxcryptset-logger")->info(
        std::string("</") + name + std::string(">\n"));

    // add it to cryptset
    //_wxcs->add(name, _tokens.back());
    _tokens.pop_back();
  }
  else if (_sym == WX_LEFT_ANGLE)
  {
    WxCryptSet _wxcs_c;
    // parse it
    _WxParser wxl(this, &_wxcs_c);
    wxl.wx_decl();
    // add it to current cryptset
    _wxcs->addSet(_wxcs_c);
    // update values of _sym and _tokenStr. This may be better
    // handled if we make _sym and _tokenStr into some sort of
    // per-class (rather than per-object) variables
    _sym = wxl._sym; _tokenStr = wxl._tokenStr;
  }
}

void
WxCryptSet::_WxParser::wx_decl_list_rest()
{
  if ((_sym == WX_ID) || (_sym == WX_LEFT_ANGLE))
  {
    wx_decl_value();
    wx_decl_list_rest();
  }
}

void
WxCryptSet::_WxParser::wx_decl_list()
{
  wx_decl_value();
  wx_decl_list_rest();
}

void
WxCryptSet::_WxParser::wx_decl()
{
  // tag open
  wx_decl_top();

  // get hold of name and set it
  _wxcs->_setName( (_tokens.back()).to_value<std::string>() );
  _tokens.pop_back();

  // inside stuff
  wx_decl_list();

  // tag close
  wx_decl_end();
}

void
WxCryptSet::_WxParser::input_file()
{
  wx_decl();
}

void
WxCryptSet::_parse(_WxParser& fl)
{
  // parse input file
  fl.input_file();
}
