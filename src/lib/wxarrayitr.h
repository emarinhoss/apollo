#ifndef __wxarrayitr__
#define __wxarrayitr__

// WarpX includes

// std includes

// forward declare the WxArrayBase class
template <typename T, int TYPE> class WxArrayBase;

/**
 * Iterator over components of data stored at an array index
 */
template<class T>
class WxArrayItr 
{
  public:
    // friend so only arraybase can instantiate
    template <typename DT, int TYPE> friend class WxArrayBase;

    WxArrayItr() {}

/**
 * Copy the iterator: a copy of the iterator points to the same
 * underlying data array
 */
    WxArrayItr(const WxArrayItr<T>& f)
            : _numComponents(f._numComponents),
              _data(f._data) {
    }

/**
 * Assign the iterator: a copy of the iterator points to the same
 * underlying data array
 */
    WxArrayItr<T>& operator=(const WxArrayItr<T>& f) {
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

/**
 * Return component of iterator
 *
 * @param k component required
 * @return value of component k
 */
    T operator[](unsigned k) const {
        return _data[k];
    }

/**
 * Return component of iterator
 *
 * @param k component required
 * @return reference to component k
 */
    T& operator[](unsigned k) {
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
 * @param numComponents number of components in array
 * @param data pointer to data iterated by iterator
 */
    WxArrayItr(unsigned numComponents, T* data)
            : _numComponents(numComponents) ,
              _data(data) {
    }

    unsigned _numComponents;
    T *_data;
};

#endif //  __wxarrayitr__
