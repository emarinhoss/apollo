#ifndef __wxboxneighholder__
#define __wxboxneighholder__

// WarpX includes
#include "wxbox.h"
#include "wxboxcmp.h"
#include "wxboxneigh.h"

// std includes
#include <map>
#include <vector>

/**
 * WxBoxNeighHolder holds a set of objects of type WxBoxNeigh each of
 * which can compute and store neighbor information for a given pad
 * cell distribution. The neighbor information is computed only once
 * for every unique pad cell distribution.
 */
template<typename TYPE>
class WxBoxNeighHolder
{
  public:
/**
 * Return the neighbox list for specified pad cell distribution
 *
 * @param n box number whose neighbor list is required
 * @param box split box whose neighbor calculation is needed
 * @param lp distribution of pad cells along lower edge in each direction
 * @param up distribution of pad cells along upper edge in each direction
 * @return list of neighboring box numbers
 */
    std::vector<unsigned> 
    getNeigh(unsigned n, const WxSplitBox<TYPE>& box, TYPE lp[], TYPE up[]);

  private:
    typedef std::map<WxBox<TYPE>,  WxBoxNeigh<TYPE>*, WxBoxCmp<TYPE> > NeighPtrMap_t;
    typedef std::pair<WxBox<TYPE>, WxBoxNeigh<TYPE>* > NeighPtrPair_t;

    NeighPtrMap_t _neighPtrMap;
};

#endif //  __wxboxneighholder__
