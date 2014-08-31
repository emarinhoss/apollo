#ifndef __wxarraybase__h__
#define __wxarraybase__h__

// WarpX includes
#include "wxarraydata.h"
#include "wxarrayitr.h"
#include "wxconstarrayitr.h"
#include "wxindexer.h"
#include "wxrange.h"
#include "wxobject.h"
#include "wxmath.h"

template <typename T, int TYPE>
class WxArrayBase
{
  public:
/** Setup the object */
    virtual void setup(const WxCryptSet& wxc);

/** Delete object */
    virtual ~WxArrayBase() {
    }

/**
 * Number of elements at each index location
 */
    unsigned numComponents() const {
      return _numComponents;
    }

/** Return x offsets for the components */
    std::vector<double> xoffsets() const {
      return _xoffsets;
    }

/** Return y offsets for the components */
    std::vector<double> yoffsets() const {
      return _yoffsets;
    }

/** Return z offsets for the components */
    std::vector<double> zoffsets() const {
      return _zoffsets;
    }

/**
 * Create a new iterator for iterating over the components stored at a
 * array location. The iterator is attached to the components at the
 * first index location in the array.
 *
 * @return new iterator object
 */
    WxArrayItr<T> createItr();

/**
 * Create a new const iterator for iterating over the components
 * stored at a array location. The iterator is attached to the
 * components at the first index location in the array.
 *
 * @return new iterator object
 */
    WxConstArrayItr<T> createConstItr() const;

/**
 * Attach iterator to the components at specified index location in
 * the array.
 *
 * @param itr iterator
 * @param indices index location to attach to
 */
    void setItr(WxArrayItr<T>& itr, int *indices);

/**
 * Attach iterator to the components at specified index location in
 * the array.
 *
 * @param itr iterator
 * @param k0 index location in direction 0
 */
    void setItr(WxArrayItr<T>& itr, int k0);

/**
 * Attach iterator to the components at specified index location in
 * the array.
 *
 * @param itr iterator
 * @param k0 index location in direction 0
 * @param k1 index location in direction 1
 */
    void setItr(WxArrayItr<T>& itr, int k0, int k1);

/**
 * Attach iterator to the components at specified index location in
 * the array.
 *
 * @param itr iterator
 * @param k0 index location in direction 0
 * @param k1 index location in direction 1
 * @param k2 index location in direction 2
 */
    void setItr(WxArrayItr<T>& itr, int k0, int k1, int k2);

/**
 * Attach iterator to the components at specified index location in
 * the array.
 *
 * @param itr iterator
 * @param k0 index location in direction 0
 * @param k1 index location in direction 1
 * @param k2 index location in direction 2
 * @param k3 index location in direction 3
 */
    void setItr(WxArrayItr<T>& itr, int k0, int k1, int k2, int k3);

/**
 * Attach iterator to the components at specified index location in
 * the array.
 *
 * @param itr iterator
 * @param indices index location to attach to
 */
    void setItr(WxConstArrayItr<T>& itr, int *indices) const;

/**
 * Attach iterator to the components at specified index location in
 * the array.
 *
 * @param itr iterator
 * @param k0 index location in direction 0
 */
    void setItr(WxConstArrayItr<T>& itr, int k0) const;

/**
 * Attach iterator to the components at specified index location in
 * the array.
 *
 * @param itr iterator
 * @param k0 index location in direction 0
 * @param k1 index location in direction 1
 */
    void setItr(WxConstArrayItr<T>& itr, int k0, int k1) const;

/**
 * Attach iterator to the components at specified index location in
 * the array.
 *
 * @param itr iterator
 * @param k0 index location in direction 0
 * @param k1 index location in direction 1
 * @param k2 index location in direction 2
 */
    void setItr(WxConstArrayItr<T>& itr, int k0, int k1, int k2) const;

/**
 * Attach iterator to the components at specified index location in
 * the array.
 *
 * @param itr iterator
 * @param k0 index location in direction 0
 * @param k1 index location in direction 1
 * @param k2 index location in direction 2
 * @param k3 index location in direction 3
 */
    void setItr(WxConstArrayItr<T>& itr, int k0, int k1, int k2, int k3) const;

/**
 * General indexing routines.
 *
 * @param indices 'indices[k]' is index into kth dimension
 * @return value at index location
 */
    T operator()(const int *indices) const;
    T& operator()(const int *indices);

/**
 * Rank-1 indexer
 * 
 * @param k0 index into dimension 0
 */
    T  operator()(int k0) const;
    T& operator()(int k0);

/**
 * Rank-2 indexer
 * 
 * @param k0 index into dimension 0
 * @param k1 index into dimension 1
 */
    T  operator()(int k0, int k1) const;
    T& operator()(int k0, int k1);

/**
 * Rank-3 indexer
 * 
 * @param k0 index into dimension 0
 * @param k1 index into dimension 1
 * @param k2 index into dimension 2
 */
    T  operator()(int k0, int k1, int k2) const;
    T& operator()(int k0, int k1, int k2);

/**
 * Rank-4 indexer
 * 
 * @param k0 index into dimension 0
 * @param k1 index into dimension 1
 * @param k2 index into dimension 2
 * @param k3 index into dimension 3
 */
    T  operator()(int k0, int k1, int k2, int k3) const;
    T& operator()(int k0, int k1, int k2, int k3);

/**
 * Pointer to underlying data: modifying this pointer directly can be
 * rather dangerous
 */
    T* data();

    unsigned _numComponents;
    WxIndexer<TYPE> _indexer;
    WxArrayData<T> _data;

  protected:
/**
 * Default ctor: do not use this directly
 */
    WxArrayBase() {
    }

/**
 * Array spanning given range object
 *
 * @param range range spanned by array
 * @param value value to assign to each element
 */
    WxArrayBase(unsigned numComponents, const WxRange& range, const T& value);

/**
 * Array spanning given range object and attached to supplied memory
 *
 * @param range range spanned by array
 * @param cdata raw data pointer for storing array data
 */
    WxArrayBase(unsigned numComponents, const WxRange& range, T *cdata);

/**
 * Copy constructor does not allocate new memory. The created array
 * shares the data with the original one.
 *
 * @param f array to copy
 */
    WxArrayBase(const WxArrayBase<T, TYPE>& f);

/**
 * Assignment operator does not allocate new memory. The created array
 * shares the data with the original one.
 *
 * @param f array to copy
 * @return reference to new array 
 */
    WxArrayBase<T, TYPE>& operator=(const WxArrayBase<T, TYPE>& f);

  private:
/* Offsets in each cell for components */
    std::vector<double> _xoffsets, _yoffsets, _zoffsets;
};

template <typename T, int TYPE>
WxArrayBase<T, TYPE>::
WxArrayBase(unsigned numComponents, const WxRange& range, const T& value)
  : _numComponents(numComponents),
    _indexer(range), 
    _data(numComponents*range.size()) 
{

  // initialize to given value
  for (unsigned i=0; i<_numComponents*range.size(); ++i)
    _data.data[i] = value;

  // set all offsets to 0.5 by default
  for (unsigned i=0; i<_numComponents; ++i)
  {
    _xoffsets.push_back(0.5);
    _yoffsets.push_back(0.5);
    _zoffsets.push_back(0.5);
  }
}

template <typename T, int TYPE>
WxArrayBase<T, TYPE>::
WxArrayBase(unsigned numComponents, const WxRange& range, T *cdata)
  : _numComponents(numComponents),
    _indexer(range),
    _data(numComponents*range.size(), cdata) 
{

  // set all offsets to 0.5 by default
  for (unsigned i=0; i<_numComponents; ++i)
  {
    _xoffsets.push_back(0.5);
    _yoffsets.push_back(0.5);
    _zoffsets.push_back(0.5);
  }
}

template <typename T, int TYPE>
WxArrayBase<T, TYPE>::
WxArrayBase(const WxArrayBase<T, TYPE>& f)
  : _numComponents(f._numComponents),
    _indexer(f._indexer), 
    _data(f._data),
    _xoffsets(f._xoffsets),
    _yoffsets(f._yoffsets),
    _zoffsets(f._zoffsets)
{
}

template <typename T, int TYPE>
WxArrayBase<T, TYPE>&
WxArrayBase<T, TYPE>::operator=(const WxArrayBase<T, TYPE>& f) {
  if (this==&f) return *this;
  _numComponents = f._numComponents;
  _indexer = f._indexer;
  _data = f._data;
  _xoffsets = f._xoffsets;
  _yoffsets = f._yoffsets;
  _zoffsets = f._zoffsets;
  return *this;
}

template <typename T, int TYPE>
void
WxArrayBase<T, TYPE>::setup(const WxCryptSet& wxc) 
{
  // read offsets for each component
  if (wxc.has("xoffsets"))
  {
    std::vector<WxAny> offsets = wxc.template get<std::vector<WxAny> >
      ("xoffsets");
    for (unsigned i=0; i< dmin<unsigned>(offsets.size(), _numComponents); ++i)
      _xoffsets[i] = wx_any_cast<double>(offsets[i]);
  }
  if (wxc.has("yoffsets"))
  {
    std::vector<WxAny> offsets = wxc.template get<std::vector<WxAny> >
      ("yoffsets");
    for (unsigned i=0; i< dmin<unsigned>(offsets.size(), _numComponents); ++i)
      _yoffsets[i] = wx_any_cast<double>(offsets[i]);
  }
  if (wxc.has("zoffsets"))
  {
    std::vector<WxAny> offsets = wxc.template get<std::vector<WxAny> >
      ("zoffsets");
    for (unsigned i=0; i< dmin<unsigned>(offsets.size(), _numComponents); ++i)
      _zoffsets[i] = wx_any_cast<double>(offsets[i]);
  }
}

template <typename T, int TYPE>
inline
WxArrayItr<T>
WxArrayBase<T, TYPE>::createItr() 
{
  return WxArrayItr<T>(_numComponents, &_data.data[0]);
}

template <typename T, int TYPE>
inline
WxConstArrayItr<T>
WxArrayBase<T, TYPE>::createConstItr() const
{
  return WxConstArrayItr<T>(_numComponents, &_data.data[0]);
}

template <typename T, int TYPE>
inline
void 
WxArrayBase<T, TYPE>::setItr(WxArrayItr<T>& itr, int *indices)
{
  itr._data = &_data.data[
      _numComponents*_indexer.index(indices)];
}

template <typename T, int TYPE>
inline
void 
WxArrayBase<T, TYPE>::setItr(WxArrayItr<T>& itr, int k0)
{
  itr._data = &_data.data[
      _numComponents*_indexer.index(k0)];
}

template <typename T, int TYPE>
inline
void 
WxArrayBase<T, TYPE>::setItr(WxArrayItr<T>& itr, int k0, int k1)
{
  itr._data = &_data.data[
      _numComponents*_indexer.index(k0,k1)];
}

template <typename T, int TYPE>
inline
void 
WxArrayBase<T, TYPE>::setItr(WxArrayItr<T>& itr, int k0, int k1, int k2)
{
  itr._data = &_data.data[
      _numComponents*_indexer.index(k0,k1,k2)];
}

template <typename T, int TYPE>
inline
void 
WxArrayBase<T, TYPE>::setItr(WxArrayItr<T>& itr, int k0, int k1, int k2, int k3)
{
  itr._data = &_data.data[
      _numComponents*_indexer.index(k0,k1,k2,k3)];
}

template <typename T, int TYPE>
inline
void
WxArrayBase<T, TYPE>::setItr(WxConstArrayItr<T>& itr, int *indices) const
{
  itr._data = &_data.data[
      _numComponents*_indexer.index(indices)];
}

template <typename T, int TYPE>
inline
void 
WxArrayBase<T, TYPE>::setItr(WxConstArrayItr<T>& itr, int k0) const
{
  itr._data = &_data.data[
      _numComponents*_indexer.index(k0)];
}

template <typename T, int TYPE>
inline
void 
WxArrayBase<T, TYPE>::setItr(WxConstArrayItr<T>& itr, int k0, int k1) const
{
  itr._data = &_data.data[
      _numComponents*_indexer.index(k0,k1)];
}

template <typename T, int TYPE>
inline
void 
WxArrayBase<T, TYPE>::setItr(WxConstArrayItr<T>& itr, int k0, int k1, int k2) const
{
  itr._data = &_data.data[
      _numComponents*_indexer.index(k0,k1,k2)];
}

template <typename T, int TYPE>
inline
void 
WxArrayBase<T, TYPE>::setItr(WxConstArrayItr<T>& itr, int k0, int k1, int k2, int k3) const
{
  itr._data = &_data.data[
      _numComponents*_indexer.index(k0,k1,k2,k3)];
}

template <typename T, int TYPE>
inline
T
WxArrayBase<T, TYPE>::operator()(const int *indices) const
{
  return _data.data[
      _numComponents*_indexer.index(indices)];
}

template <typename T, int TYPE>
inline
T& 
WxArrayBase<T, TYPE>::operator()(const int *indices)
{
  return _data.data[
      _numComponents*_indexer.index(indices)];
}

template <typename T, int TYPE>
inline
T
WxArrayBase<T, TYPE>::operator()(int k1) const
{
  return _data.data[
      _numComponents*_indexer.index(k1)]; 
}

template <typename T, int TYPE>
inline
T&
WxArrayBase<T, TYPE>::operator()(int k1)
{ 
  return _data.data[
      _numComponents*_indexer.index(k1)];
}

template <typename T, int TYPE>
inline
T 
WxArrayBase<T, TYPE>::operator()(int k1, int k2) const
{
  return _data.data[
      _numComponents*_indexer.index(k1,k2)];
}

template <typename T, int TYPE>
inline
T& 
WxArrayBase<T, TYPE>::operator()(int k1, int k2)
{
  return _data.data[
      _numComponents*_indexer.index(k1,k2)]; 
}

template <typename T, int TYPE>
inline
T 
WxArrayBase<T, TYPE>::operator()(int k1, int k2, int k3) const
{
  return _data.data[
      _numComponents*_indexer.index(k1,k2,k3)]; 
}

template <typename T, int TYPE>
inline
T& 
WxArrayBase<T, TYPE>::operator()(int k1, int k2, int k3)
{
  return _data.data[
      _numComponents*_indexer.index(k1,k2,k3)]; 
}

template <typename T, int TYPE>
inline
T 
WxArrayBase<T, TYPE>::operator()(int k1, int k2, int k3, int k4) const
{
  return _data.data[
      _numComponents*_indexer.index(k1,k2,k3,k4)]; 
}

template <typename T, int TYPE>
inline
T& 
WxArrayBase<T, TYPE>::operator()(int k1, int k2, int k3, int k4)
{
  return _data.data[
      _numComponents*_indexer.index(k1,k2,k3,k4)]; 
}

template <typename T, int TYPE>
inline
T*
WxArrayBase<T, TYPE>::data() 
{ 
  return this->_data.data; 
}

#endif //  __wxarraybase__h__
