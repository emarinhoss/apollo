// WarpX includes
#include "wxcryptsetlexer.h"
#include "wxexcept.h"

// std includes
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <cstdlib>
#include <cstdio>

#include <string.h>
#include <ctype.h>

int
WxCryptSetLexer::_yylex()
{
  char c;
  // skip white spaces, new lines  and tabs
  while ((c=_is.get()) == ' ' || c == '\t')
    ;

  // end-of-file
  if (c == EOF)
    return WX_DONE;
  // line continuation
  if (c == '\\')
  {
    c = _is.get();
    if (c == '\n')
    {
      _lineno++;
      return yylex();
    }
  }
  // comment
  if (c == '#')
  {
    while ((c=_is.get()) != '\n' && c != EOF)
      ;
    if (c=='\n')
      _lineno++;
    return _yylex();
  }
  // number
  if (c == '-' || c == '+' || c == '.' || isdigit(c))
  {
    _yytext = ""; // clear token string
    _is.putback(c);
    return _scan_number();
  }
  // name
  if (isalpha(c) || c == '_' || c == '/' || c == '~')
  {
    _yytext = ""; // clear token string
    do
    {
      _yytext.push_back(c);
    } while ((c=_is.get()) != EOF && (_isidchar(c)));
    _is.putback(c);
    return WX_KEY;
  }
  // string
  if (c == '"') 
  {
    _yytext = ""; // clear token string
    while ((c=_is.get()) != '"')
    {
      if (c == '\n' || c == EOF) {
        WxExcept wxe;
        wxe << "Missing quote for string";
        throw wxe;
      }
      _yytext.push_back(this->_backslash(c));
    }
    return WX_STRING;
  }
  // operators, assignment, definitions and declaration
  switch (c) 
  {
    case '>':	return WX_RIGHT_ANGLE;
    case '<':	return _follow('/', WX_LEFT_ANGLE_FRONT_SLASH, WX_LEFT_ANGLE);
    case '=':	return WX_EQUAL;
    case ',':	return WX_COMMA;
    case '[':	return WX_LEFT_BOX;
    case ']':	return WX_RIGHT_BOX;
    case ';':       return WX_SEMI_COLON;
    case '\n':	_lineno++; return _yylex(); // skip over a newline
    default:	return c;
  }
}

char 
WxCryptSetLexer::_backslash(char c) 
{
  static char transtab[] = "b\bf\fn\nr\rt\t";
  if (c != '\\')
    return c;
  c = _is.get();
  if (islower(c) && strchr(transtab, c))
    return strchr(transtab, c)[1];
  return c;
}

unsigned 
WxCryptSetLexer::_follow(char expect, unsigned ifyes, unsigned ifno) 
{
  int c = _is.get();
  // check next character character
  if (c == expect) return ifyes;
  // unexpected one
  _is.putback(c);
  return ifno;
}

bool
WxCryptSetLexer::_isidchar(char c)
{
  return (isalnum(c) ||
    c == '_' || 
    c == '-' ||
    c == '.' ||
    c == '/' ||
    c == '~' ||
    c == ':');
}

// NOTE: this function needs much improvement
int
WxCryptSetLexer::_scan_number()
{
  int state = 0; // start state
  bool isint = true; // true of number scanned is an integer
  char c;

  // following only recognizes if the token is a number. the actual
  // conversion is handled by atof and atoi (FIX THIS: MAIN ISSUE IS
  // THAT AN OVERFLOW MAY OCCUR AND HANDLING THIS IS TRICKY)
  while( 1 )
  {
    c = _is.get();
    switch(state)
    {
      case 0:
        // +,-,.,[0-9]
        if (c == '+')
          // ignore +ve sign
          state = 1;
        else if(c == '-')
          state = 1;
        else if (isdigit(c))
          state = 1;
        else if(c == '.')
          state = 2;
        else
          // error
          state = -1;

        break;

      case 1:
        // [0-9], ., [e,E]
        if (isdigit(c))
          state = 1;
        else if (c == '.')
          state = 2;
        else if((c == 'e') || (c == 'E'))
          state = 3;
        else
          // error
          state = -1;

        break;

      case 2:
        isint = false;
        // [0-9], [e,E]
        if (isdigit(c))
          state = 1;
        else if ((c == 'e') || (c == 'E'))
          state = 3;
        else
          state = -1;

        break;

      case 3:
        isint = false;
        // [+,-], [0-9]
        if(isdigit(c))
          state = 4;
        else if( c == '+')
          state = 4;
        else if (c == '-')
          state = 4;
        else
          state = -1;

        break;

      case 4:
        isint = false;
        // [0-9]
        if(isdigit(c))
          state = 4;
        else
          state = -1;
    }
    if(state == -1)
      break; // finish parsing
    _yytext.push_back(c); // push character to current token string
  }
  _is.putback(c);
  // compute number
  if (isint)
  {
    _integer = atoi(_yytext.c_str());
    return WX_INT;
  }
  else
  {
    _real = atof(_yytext.c_str());
    return WX_REAL;
  }
}
