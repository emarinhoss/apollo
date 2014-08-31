/**
 * @file	wxexcept.cc
 *
 * @brief	Class to represent exceptions in Warpx
 *
 * @version	$Id: wxexcept.cc 4080 2009-02-08 23:43:19Z AmmarHakim $
 *
 * Copyright &copy; 2008-2009, Ammar Hakim. Released under Eclipse
 * Licence version 1.0.
 */

// lib includes
#include <wxexcept.h>

namespace Warpx
{
  Except::Except() 
  {
  }

  Except::Except(const std::string& str)
    : std::exception() 
  {
    (*this) << str; 
  }

  Except::Except(const Except& ex)
  {
    (*this) << ex.exceptStrm.str();
  }

  Except::~Except() throw()
  {
  }

  Except&
  Except::operator=(const Except& ex)
  {
    if (this==&ex) return *this;
    this->exceptStrm.str(""); // zap the string
    (*this) << ex.exceptStrm.str();
    return *this;
  }

  const char*
  Except::what() const throw()
  {
    return this->exceptStrm.str().c_str();
  }
}
