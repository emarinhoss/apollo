#ifndef __wxcryptsetlexer__h__
#define __wxcryptsetlexer__h__

// WarpX includes

// std includes
#include <iostream>
#include <istream>
#include <string>

// WarpX input file lexer tokens
enum {
  WX_ERROR = -2,
  WX_DONE = -1,
  WX_LEFT_BOX = 256,
  WX_RIGHT_BOX,
  WX_LEFT_ANGLE,
  WX_RIGHT_ANGLE,
  WX_ID,
  WX_KEY = WX_ID, // key and value are just IDs
  WX_VALUE = WX_ID,
  WX_EQUAL,
  WX_FRONT_SLASH,
  WX_COMMA,
  WX_INT,
  WX_REAL,
  WX_STRING,
  WX_SEMI_COLON,
  WX_LEFT_ANGLE_FRONT_SLASH
};

class WxCryptSetLexer
{
  public:
/**
 * Set input stream to one supplied. This defaults to standard input.
 *
 * @param is Input stream to scan for characters
 */
    WxCryptSetLexer(std::istream& is)
      : _is(is), _lineno(1), _integer(0), _real(0.0) {
    }

/**
 * Scans input stream and returns a single token.
 */
    int yylex() {
      _lastSym = _yylex();
      return _lastSym;
    }

/**
 * Returns the current line number being scanned
 */
    unsigned lineno() const {
      return _lineno;
    }

/**
 * Returns a character pointer representing the current token scanned.
 */
    std::string YYText() const {
      return _yytext;
    }

/**
 * Returns integer scanned
 */
    int integer() const {
      return _integer;
    }

/**
 * Returns real number scanned
 */
    double real() const {
      return _real;
    }

    int lastSym() const {
      return _lastSym;
    }

  private:
    std::istream& _is; // input stream to tokenize
    unsigned _lineno;
    std::string _yytext; // text of last token read
    int _integer; // scanned integer
    double _real; // scanned real
    unsigned _lastSym; // last symbol scanned

    // do the actual tokenization
    int _yylex();

    // get next char with \'s interpreted
    char _backslash(char c);

    // look ahead for two character operators
    unsigned _follow(char expect, unsigned ifyes, unsigned ifno);

    // check if character can belong to WX_ID
    bool _isidchar(char c);

    // scans a number
    int _scan_number();
};

#endif //  __wxcryptsetlexer_h__
