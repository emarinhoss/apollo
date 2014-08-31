#include "wxcrypt.h"


WxCrypt::WxCrypt() 
{
}

WxCrypt::~WxCrypt() 
{
  // delete all entries in map
  _values.erase( _values.begin(), _values.end() );
}

WxCrypt::WxCrypt(const WxCrypt& crypt)
{
  AnyMap_t::const_iterator i;
  // copy all entries
  for (i=crypt._values.begin(); i!=crypt._values.end(); ++i)
    _values.insert( AnyPair_t((*i).first, (*i).second) );
  // copy the type -> names map
  copyNames<int>(crypt);
  copyNames<double>(crypt);
  copyNames<std::string>(crypt);
  copyNames<WxAny>(crypt);
  copyNames<std::vector<WxAny> >(crypt);
}

WxCrypt&
WxCrypt::operator=(const WxCrypt& rhs)
{
  if (this==&rhs) return *this;

  // delete all entries in map
  _values.erase( _values.begin(), _values.end() );
  // delete all entries in type -> name map
  wxTypeMapExtract<int>(_typeToNames).names.erase(
      wxTypeMapExtract<int>(_typeToNames).names.begin(),
      wxTypeMapExtract<int>(_typeToNames).names.end());

  wxTypeMapExtract<int>(_typeToNames).names.erase(
      wxTypeMapExtract<int>(_typeToNames).names.begin(),
      wxTypeMapExtract<int>(_typeToNames).names.end());

  wxTypeMapExtract<double>(_typeToNames).names.erase(
      wxTypeMapExtract<double>(_typeToNames).names.begin(),
      wxTypeMapExtract<double>(_typeToNames).names.end());

  wxTypeMapExtract<std::string>(_typeToNames).names.erase(
      wxTypeMapExtract<std::string>(_typeToNames).names.begin(),
      wxTypeMapExtract<std::string>(_typeToNames).names.end());

  wxTypeMapExtract<WxAny>(_typeToNames).names.erase(
      wxTypeMapExtract<WxAny>(_typeToNames).names.begin(),
      wxTypeMapExtract<WxAny>(_typeToNames).names.end());

  wxTypeMapExtract<std::vector<WxAny> >(_typeToNames).names.erase(
      wxTypeMapExtract<std::vector<WxAny> >(_typeToNames).names.begin(),
      wxTypeMapExtract<std::vector<WxAny> >(_typeToNames).names.end());

  // add entries from rhs
  AnyMap_t::const_iterator i;
  for (i=rhs._values.begin(); i!=rhs._values.end(); ++i)
    _values.insert( AnyPair_t((*i).first, (*i).second) );

  // copy the type -> names map
  copyNames<int>(rhs);
  copyNames<double>(rhs);
  copyNames<std::string>(rhs);
  copyNames<std::vector<WxAny> >(rhs);

  return *this;
}

bool 
WxCrypt::has(const std::string& name) const
{
  AnyMap_t::const_iterator i;
  i = _values.find(name);
  return (i != _values.end()) ? true : false;
}
