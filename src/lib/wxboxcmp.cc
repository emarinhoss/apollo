#include "wxboxcmp.h"

template <typename TYPE>
bool
WxBoxCmp<TYPE>::operator()(const WxBox<TYPE>& lhs, const WxBox<TYPE>& rhs) const
{
    for (unsigned i=0; i<lhs.ndims(); ++i)
    {
        if (lhs.lower(i) < rhs.lower(i))
            return true;
        else if (lhs.lower(i) > rhs.lower(i))
            return false;
    }
    for (unsigned i=0; i<lhs.ndims(); ++i)
    {
        if (lhs.upper(i) < rhs.upper(i))
            return true;
        else if (lhs.upper(i) > rhs.upper(i))
            return false;
    }
    return false; // lhs is identical to rhs
}

// instantiations
template class WxBoxCmp<int>;
//template class WxBoxCmp<float>;
template class WxBoxCmp<double>;
