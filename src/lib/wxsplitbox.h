#ifndef __wxsplitbox__
#define __wxsplitbox__

// WarpX lib includes
#include "wxbox.h"
#include "wxboxneighholder.h"
#include "wxmsgbase.h"

// std includes
#include <map>
#include <vector>

/**
 * WxSplitBox splits a box into smaller boxes
 */
template<typename TYPE>
class WxSplitBox
{
  public:

/**
 * Split supplied box into smaller number of boxes. Each sub-box is
 * labelled from 0..nboxes-1
 *
 * @param msgBase messenger associated with this decomposition
 * @param box box to split
 * @param pdirs pdisr[i] is true if domain is periodic in direction 'i'
 */
    WxSplitBox(WxMsgBase *msgBase, const WxBox<TYPE>& box, bool pdirs[]=0);

/**
 * Split supplied box into smaller boxes. Each sub-box is labelled
 * from [0..nbox), where nbox is the number of sub-boxes. The number
 * of boxes must be less than the number of processes in the messenger
 * returned by msgBase->numProcs().
 *
 * @param msgBase messenger associated with this decomposition
 * @param subBox[] number of sub-boxes along dimension dim is subBox[dim]
 * @param box global box to split
 * @param pdirs pdisr[i] is true if domain is periodic in direction 'i'
 */
    WxSplitBox(WxMsgBase *msgBase, unsigned subBox[], const WxBox<TYPE>& box, bool pdirs[]=0);

/**
 * Destory a split box
 */
    virtual ~WxSplitBox();

/**
 * Copy ctor makes a shallow copy of the supplied box
 */
    WxSplitBox(const WxSplitBox<TYPE>& bs);

/**
 * Assignment operator makes a shallow copy of the supplied box
 */
    WxSplitBox<TYPE>& operator=(const WxSplitBox<TYPE>& bs);

/**
 * Dimension of box split
 */
    unsigned ndims() const {
      return _data->box.ndims();
    }

/**
 * Number of sub-boxes in split
 */
    unsigned nboxes() const {
      return _data->nboxes;
    }

/**
 * Number of sub-boxes including pseudo halo-boxes (used only in cases
 * of periodic boundaries)
 *
 * @return number of sub-boxes in the system
 */
    unsigned getNumPseudoRanks() const {
      return _data->boxMap.size();
    }

/**
 * Return the global box that was split.
 */
    WxBox<TYPE> getGlobalBox() const {
      return _data->box;
    }

/**
 * Return the messenger for this split.
 *
 * @return pointer to WxMsgBase object which this split box was
 * constructed with.
 */
    WxMsgBase* getMsgr() const {
      return _data->msgBase;
    }

/**
 * Return neighbors for a given pad cell distribution
 *
 * @param n box number for which neighbor information is needed
 * @param lp distribution of pad cells along lower edge in each direction
 * @param up distribution of pad cells along upper edge in each direction
 * @return list of neighboring box numbers
 */
    std::vector<unsigned> getNeigh(unsigned n, TYPE lp[], TYPE up[]) {
      return _data->neighHolder->getNeigh(n, *this, lp, up);
    }

/**
 * Return a sub-box in the split.
 *
 * @param n sub-box number. Should be in [0..nboxes)
 * @return sub-box associated with 'n'
 */
    WxBox<TYPE> getBox(unsigned n) const;

/**
 * Given a box number, return the rank on which it exists
 *
 * @param n the box number
 * @return the rank on which the box (or pseudobox) exists
 */
    unsigned getBoxRank(unsigned n) const;

/**
 * Given a point, return which box it resides in
 *
 * @param coords coordinates of the point in question
 * @return rank of the box in which the point resides
 */
    unsigned getRank(TYPE coords[]) const;

  private:

    // types of unsigned -> WxBox map and pairs
    typedef std::map<unsigned, WxBox<TYPE> > BoxMap_t;
    typedef std::pair<unsigned, WxBox<TYPE> > BoxPair_t;
    typedef std::map<unsigned, unsigned> BoxRank_t;
    typedef std::pair<unsigned, unsigned> BoxRPair_t;

    // private struct to hold split-box data
    struct WxSplitBoxData
    {
        WxSplitBoxData(unsigned nboxes, const WxBox<TYPE>& box, WxMsgBase *msgBase)
          : nboxes(nboxes),
            periods(new bool[box.ndims()]),
            box(box), 
            neighHolder(new WxBoxNeighHolder<TYPE>),
            msgBase(msgBase) {
        }

        virtual ~WxSplitBoxData() {
          delete neighHolder;
        }

        unsigned nboxes;
        bool *periods;
        WxBox<TYPE> box;
        WxBoxNeighHolder<TYPE> *neighHolder;
        WxMsgBase *msgBase;
        BoxMap_t boxMap;
        BoxRank_t boxRank;
    };

    unsigned NDIMS;
    int *_useCount;
    WxSplitBoxData *_data;

    void _divideBox(unsigned num, const WxBox<TYPE>& decompBox);
    void _duplicateBoxes();
    std::vector<WxBox<TYPE> > _breakBoxes(unsigned subBox[16], const WxBox<TYPE>& box);
};

#endif //  __wxsplitbox__
