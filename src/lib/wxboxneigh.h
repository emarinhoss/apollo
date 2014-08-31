#ifndef __wxboxneigh__
#define __wxboxneigh__

// WarpX lib includes
#include "wxbox.h"

// std includes
#include <vector>
#include <map>

// forward declare split box
template <typename T> class WxSplitBox;

template<typename TYPE>
class WxBoxNeigh
{
  public:
/**
 * Given a boxsplitter object and a distribution of pad cells creates
 * neighbor links.
 */
    WxBoxNeigh(const WxSplitBox<TYPE> &ws, TYPE lp[], TYPE up[]);

/**
 * Neighbor box numbers of the specified box number
 */
    std::vector<unsigned> getNeigh(unsigned n) const;

  private:
    typedef std::map<unsigned, std::vector<unsigned> > NeighMap_t;
    typedef std::pair<unsigned, std::vector<unsigned> > NeighPair_t;

    TYPE lp[16], up[16]; // distribution of pad cells
    NeighMap_t _neighMap; // map of a box number to number of its neighboring boxes
};

#endif //  __wxboxneigh__
