#include "wxsplitbox.h"
#include "wxboxneighholder.h"

template <typename TYPE>
std::vector<unsigned> 
WxBoxNeighHolder<TYPE>::
getNeigh(unsigned n, const WxSplitBox<TYPE>& box, TYPE lp[], TYPE up[])
{
    // first check if this pad cell distribution already occurs
    WxBox<TYPE> tmp(box.ndims(), lp, up);
    typename NeighPtrMap_t::iterator i = _neighPtrMap.find(tmp);
    if (i != _neighPtrMap.end())
        return (i->second)->getNeigh(n);

    // it does not so compute it and insert it
    WxBoxNeigh<TYPE> *ptr = new WxBoxNeigh<TYPE>(box, lp, up);
    _neighPtrMap.insert(
        NeighPtrPair_t(tmp, ptr) );

    // return newly computed neighbors
    return ptr->getNeigh(n);
}

// instantiations
// instantiations
template class WxBoxNeighHolder<int>;
//template class WxBoxNeighHolder<float>;
template class WxBoxNeighHolder<double>;

