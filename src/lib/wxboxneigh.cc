// WarpX lib includes
#include "wxboxneigh.h"
#include "wxsplitbox.h"
#include "wxexcept.h"

template<typename TYPE>
WxBoxNeigh<TYPE>::WxBoxNeigh(const WxSplitBox<TYPE> &ws, TYPE lp[], TYPE up[])
{
    // a box's neighbors are defined as those boxes which intersect
    // its extended region. STRICTLT SPEAKING THESE ARE RECIEVE
    // NEIGHBORS, I.E. GUYS FROM WHICH WE RECIEVE DATA. IN WEIRD
    // SENARIOS THESE MAY NOT BE THE SAME AS THE ONES WE SEND DATA
    // TO. WORRY ABOUT THIS LATER
    for (unsigned i=0; i<ws.nboxes(); ++i)
    {
        WxBox<TYPE> myBox = ws.getBox(i).extend(lp, up);
        std::vector<unsigned> neigh;
        for (unsigned j=0; j<ws.getNumPseudoRanks(); ++j)
        {
            if (i==j) continue;
            WxBox<TYPE> otherBox = ws.getBox(j);
            if ( !myBox.intersect(otherBox).isEmpty() )
                neigh.push_back(j);
        }
        _neighMap.insert(
            NeighPair_t(i, neigh) );
    }
}

template<typename TYPE>
std::vector<unsigned>
WxBoxNeigh<TYPE>::getNeigh(unsigned n) const
{
    typename NeighMap_t::const_iterator i = _neighMap.find(n);
    if (i != _neighMap.end())
        return (*i).second;
    // not found: throw a runtime_error. Maybe do something better?
    WxExcept wxe;
    wxe << "Box not found";
    throw wxe;
}

// instantiations
template class WxBoxNeigh<int>;
//template class WxBoxNeigh<float>;
template class WxBoxNeigh<double>;
