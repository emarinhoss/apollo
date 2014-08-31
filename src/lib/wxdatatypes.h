#ifndef __wxdatatypes__
#define __wxdatatypes__

// WarpX includes
#include "wxtypelist.h"
#include "wxany.h"

// std includes
#include <string>
#include <vector>

// typelist for supported I/O and message-ing types. These can be
// augment if needed.
//
// NOTE: If adding more types change typelist length. Also ensure no
// duplicates exist in the list.

typedef WX_TYPELIST_18(
    bool,
    char,
    unsigned char,
    short,
    unsigned short,
    int,
    unsigned int,
    long,
    unsigned long,
    float,
    double,
    long double,
    long long int,
    WxAny,
    std::vector<WxAny>,
    std::string,
    float*,
    double*
                       ) WxDataTypes_t;

#endif //  __wxdatatypes__
