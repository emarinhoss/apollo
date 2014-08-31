/* A Flex lexical analyser specification for WarpX input files.
 */

%{
// WarpX includes
#include "wxlextok.h"

// std includes
    
%}

%option yylineno
%option noyywrap
%option c++
%s EQUAL

%%

[ \t]
#.* 

\<  return WX_LEFT_ANGLE;

\<\/ return WX_LEFT_ANGLE_FRONT_SLASH;

\>  return WX_RIGHT_ANGLE;

=   {BEGIN EQUAL; return WX_EQUAL;}

\[  return WX_LEFT_BOX;

\]  return WX_RIGHT_BOX;

,   return WX_COMMA;

;   return WX_SEMI_COLON;

[a-zA-Z_][a-zA-Z0-9\-_]* return WX_ID;

[-+]?[0-9]+ return WX_INT;

[-+]?[0-9]*([eE]?[-+]?[0-9]+) |
[-+]?[0-9]*\.[0-9]*([eE]?[-+]?[0-9]+)? return WX_REAL;

\"[^"\n]*[\"\n] return WX_STRING;

<EQUAL>[a-zA-Z_\.\/:][a-zA-Z0-9\-_\\\/\.:]* {return WX_VALUE; }
<EQUAL>\n BEGIN 0;

[ \t\n]

%%
