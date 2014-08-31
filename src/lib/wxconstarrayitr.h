#ifndef __wxconstarrayitr__
#define __wxconstarrayitr__

// WarpX includes

// std includes

// forward declare the WxArrayBase class
template <typename T, int TYPE> class WxArrayBase;

/**
 * Iterator over components of data stored at an constarray index
 */
template<class T>
class WxConstArrayItr 
{
  public:
    // friend so only constarraybase can instantiate
    template <typename DT, int TYPE> friend class WxArrayBase;

/**
 * Copy the iterator: a copy of the iterator points to the same
 * underlying data constarray
 */
    WxConstArrayItr(const WxConstArrayItr<T>& f)
      : _numComponents(f._numComponents),
        _data(f._data) {
    }

/**
 * Assign the iterator: a copy of the iterator points to the same
 * underlying data constarray
 */
    WxConstArrayItr<T>& operator=(const WxConstArrayItr<T>& f) {
      if (this==&f) return *this;
      _numComponents = f._numComponents;
      _data = f._data;
      return *this;
    }

/**
 * Return number of components in iterator
 *
 * @return components in iterator
 */
    unsigned numComponents() const {
      return _numComponents;
    }

    T operator[](unsigned k) {
      return _data[k];
    }

/**
 * Return pointer to raw data
 */
    T* data() {
      return _data;
    }

  private:

/**
 * Create a new iterator having supplied number of components attached
 * to given data pointer.
 *
 * @param numComponents number of components in constarray
 * @param data pointer to data iterated by iterator
 */
    WxConstArrayItr(unsigned numComponents, T* data)
      : _numComponents(numComponents) ,
        _data(data) {
    }

    unsigned _numComponents;
    T *_data;
};

#endif //  __wxconstarrayitr__
