#ifndef __wxgridrange__
#define __wxgridrange__

// WarpX includes
#include "wxrange.h"

// std includes

class WxGridRange
{
  public:

/**
 * Dimension of grid
 */
    unsigned ndims() const {
      return _fr.ndims();
    }

/**
 * No. of interior cells along dimension 'dim'
 */
    unsigned ncells(unsigned dim) const {
      return _r.length(dim);
    }

/**
 * No of pad cells around the box
 */
    unsigned npad() const {
      return _npad;
    }

/**
 * Range of cells covered by grid interior
 */
    const WxRange& range() const {
      return _r;
    }

/**
 * Range of cells covered by grid, including ghost cells
 */
    const WxRange& fullRange() const {
      return _fr;
    }

/**
 * Lower index along dimension 'dim'
 */
    int lowerIdx(unsigned dim) const {
      return _r.lower(dim);
    }

/**
 * Upper index along dimension 'dim'
 */
    int upperIdx(unsigned dim) const {
      return _r.upper(dim);
    }

/**
 * Index of first pad cell along dimension 'dim'
 */
    int lowerPadIdx(unsigned dim) const {
      return _fr.lower(dim);
    }

/**
 * Index of last pad cell along dimension 'dim'
 */
    int upperPadIdx(unsigned dim) const {
      return _fr.upper(dim);
    }

/**
 * Returns range object representing grid-shell along lower face of
 * dimension 'dim'.
 */
    WxRange lowerShell(unsigned dim) const;

/**
 * Returns range object representing grid-shell along upper face of
 * dimension 'dim'.
 */
    WxRange upperShell(unsigned dim) const;

/**
 * Returns range object representing ghost cells along lower face of
 * dimension 'dim'.
 */
    WxRange lowerGhost(unsigned dim) const;

/**
 * Returns range object representing ghost cells along upper face of
 * dimension 'dim'.
 */
    WxRange upperGhost(unsigned dim) const;

/**
 * Returns range object representing a layer of cells 'n' thick along
 * lower face of dimension 'dim'.
 */
    WxRange lowerLayer(unsigned dim, unsigned n) const;

/**
 * Returns range object representing a layer of cells 'n' thick along
 * upper face of dimension 'dim'.
 */
    WxRange upperLayer(unsigned dim, unsigned n) const;

/**
 * Returns range object representing a layer of cells 'n' thick along
 * lower face of dimension 'dim'. This layer also includes the ghost
 * cells in the orthogonal directions.
 */
    WxRange lowerLayerWithGhost(unsigned dim, unsigned n) const;

/**
 * Returns range object representing a layer of cells 'n' thick along
 * upper face of dimension 'dim'. This layer also includes the ghost
 * cells in the orthogonal directions.
 */
    WxRange upperLayerWithGhost(unsigned dim, unsigned n) const;

/**
 * Returns total number of cells (interior+ghost) in grid
 */
    unsigned totalCells() const {
      return _fr.size();
    }

/**
 * Returns total number of cells (interior+ghost) in grid
 */
    unsigned totalInteriorCells() const {
      return _r.size();
    }

  protected:

    // default ctor does nothing
    WxGridRange() {
    }

/**
 * Constructs a new grid range object
 *
 * @param in Range object for interior
 * @param npad No of pad cells around each side of grid
 */
    WxGridRange(const WxRange& in, unsigned npad);

    WxGridRange(const WxGridRange& wr);
    WxGridRange& operator=(const WxGridRange& wr);

    void setRanges(unsigned npad, const WxRange& r, const WxRange& fr);

  private:
    unsigned _npad; // no of pad cells
    WxRange _r, _fr; // ranges for interior and full domain
};

#endif //  __wxgridrange__
