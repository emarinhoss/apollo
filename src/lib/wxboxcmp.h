#ifndef __wxboxcmp__
#define __wxboxcmp__

// WarpX lib includes
#include "wxbox.h"

// std includes
#include <algorithm>

/**
 * Class (functor) to compare two boxes. Note that this comparison is
 * ad-hoc and mainly used to provide a rather artificial
 * lexicograhical comparison so boxes can be used as keys in a
 * std::map
 */
template <typename TYPE>
struct WxBoxCmp : public std::binary_function<WxBox<TYPE>, WxBox<TYPE>, bool>
{
/**
 * Returns true if lhs<rhs and false otherwise
 */
    bool operator()(const WxBox<TYPE>& lhs, const WxBox<TYPE>& rhs) const;
};

#endif //  __wxboxcmp__
