#ifndef __wxsequencer__h__
#define __wxsequencer__h__

// WarpX includes
#include "wxrange.h"
#include "_wxstorageorder.h"

/**
 * WxSequencer provides a way to loop over all elements of a array
 * indexed by a given range array.
 */
template<int TYPE = _WX_COL_MAJOR_ORDER>
class WxSequencer
{
  public:

/**
 * Constructs sequencer given range of domain to sequence
 *
 * @param r Range object for domain
 */
    WxSequencer(const WxRange& r);

    ~WxSequencer() {
      delete [] _indices;
    }

/**
 * Moves sequencer by one element. Returns false when there are no
 * more elements.
 */
    bool step();

/**
 * Resets sequencer
 */
    void reset() {
      _first = true;
      for (unsigned i=0; i<_n; ++i) _indices[i] = _range.lower(i);
    }

/**
 * Returns index of current element
 */
    int* indices() const {
      return _indices;
    }

  private:
    WxRange _range;
    unsigned _n;
    bool _first; 
    int *_indices;
    bool _empty;

    // no copying allowed
    WxSequencer(const WxSequencer&);
    WxSequencer& operator=(const WxSequencer&);
};

template<int TYPE>
WxSequencer<TYPE>::WxSequencer(const WxRange& r)
  : _range(r), _n(r.ndims()), _first(true), _indices(new int[r.ndims()]), _empty(false)
{
  _empty = _range.area() == 0 ? true : false;
  for (unsigned i=0; i<_n; ++i) _indices[i] = r.lower(i);
}


/**
 * Specialized for row major sequencing. In row-major sequencing the
 * rows are kept contigously and so it is more efficient if the last
 * index moves fastest.
 */
template<>
inline
bool
WxSequencer<_WX_ROW_MAJOR_ORDER>::step()
{
  if (_empty)
  { // if box sequenced has zero volume, do nothing
    return false;
  }

  if (_first)
  { // first time around: indices already set in ctor
    _first = false;
    return true;
  }

  for (int i=_n-1; i>=0; --i)
  {
    _indices[i] += 1;
    if (_indices[i] > _range.upper(i)-1)
      _indices[i] = _range.lower(i);
    else
      return true;
  }
  return false;
}

/**
 * Specialized for column major sequencing. In column-major sequencing
 * the columns are kept contigously and so it is more efficient if the
 * first index moves fastest.
 */
template<>
inline
bool
WxSequencer<_WX_COL_MAJOR_ORDER>::step()
{
  if (_empty)
  { // if box sequenced has zero volume, do nothing
    return false;
  }

  if (_first)
  { // first time around: indices already set in ctor
    _first = false;
    return true;
  }

  for (unsigned i=0; i<_n; ++i)
  {
    _indices[i] += 1;
    if (_indices[i] > _range.upper(i)-1)
      _indices[i] = _range.lower(i);
    else
      return true;
  }
  return false;
}

#endif // __wxsequencer__h__
