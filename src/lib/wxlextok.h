#ifndef __wxlextok__h__
#define __wxlextok__h__

// Tokens for use in lexical analyzer for WarpX input files.

enum {
  WX_LEFT_BOX = 1,
  WX_RIGHT_BOX,
  WX_LEFT_ANGLE,
  WX_RIGHT_ANGLE,
  WX_KEY,
  WX_VALUE,
  WX_EQUAL,
  WX_ID,
  WX_FRONT_SLASH,
  WX_COMMA,
  WX_INT,
  WX_REAL,
  WX_STRING,
  WX_SEMI_COLON,
  WX_LEFT_ANGLE_FRONT_SLASH
};

#endif //  __wxlextok__h__
