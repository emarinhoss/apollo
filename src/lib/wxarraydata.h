#ifndef __wxarraydata__
#define __wxarraydata__

/**
 * Provides basic infrastructure for array classes in WarpX
 */

// following macros are used to manipulate WxArrayData.traits
static unsigned int ar_masks[] = 
{ 0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80 };

// is array data contiguous?
#define WX_AR_CONTIGUOUS ar_masks[0]
#define WX_AR_SET_CONTIGUOUS(bit) (bit) |= WX_AR_CONTIGUOUS
#define WX_AR_CLEAR_CONTIGUOUS(bit) (bit) &= ~WX_AR_CONTIGUOUS
#define WX_AR_IS_CONTIGUOUS(bit) (bit) & WX_AR_CONTIGUOUS

// was array data allocated?
#define WX_AR_ALLOC ar_masks[1]
#define WX_AR_SET_ALLOC(bit) (bit) |= WX_AR_ALLOC
#define WX_AR_CLEAR_ALLOC(bit) (bit) &= ~WX_AR_ALLOC
#define WX_AR_IS_ALLOC(bit) (bit) & WX_AR_ALLOC

template <typename T>
struct WxArrayData 
{
    WxArrayData()
      : nelems(0), data(0),
        useCount(new int(0)), traits(0) {
    }

/**
 * Create a new workspace for storing array data.
 *
 * @param nelems no of elements in array
 */
    WxArrayData(unsigned nelems)
      : nelems(nelems), data(new T[nelems]), 
        useCount(new int(1)), traits(0) {
      WX_AR_SET_ALLOC(traits);
      WX_AR_SET_CONTIGUOUS(traits);
    }

/**
 * Create a new workspace for storing array data. No memory is
 * allocated but is attached to the supplied data pointer.
 *
 * @param nelems no of elements in array
 * @param data buffer to use to element storage
 */
    WxArrayData(unsigned nelems, T *data)
      : nelems(nelems), data(data), 
        useCount(new int(1)), traits(0) {
      WX_AR_CLEAR_ALLOC(traits);
      WX_AR_SET_CONTIGUOUS(traits);
    }

/**
 * Delete memory if no one else is refering to it
 */
    virtual ~WxArrayData() {
      if ((--*useCount == 0) && WX_AR_IS_ALLOC(traits))
      {
        delete [] data;
        delete useCount;
      }
    }

/**
 * Construct new array data from suppiled one. This is a shallow copy
 * with pointer semantics.
 *
 * @param ad array data to copy from
 */
    WxArrayData(const WxArrayData<T>& ad) 
      : nelems(ad.nelems), data(ad.data), traits(ad.traits) {
      WX_AR_CLEAR_ALLOC(traits);
      ++*ad.useCount;
      useCount = ad.useCount;
    }

/**
 * Copy array data from suppiled one. This is a shallow copy with
 * pointer semantics.
 *
 * @param ad array data to copy from
 * @return reference to newly created array data
 */
    WxArrayData<T>& operator=(const WxArrayData<T>& ad) {
      ++*ad.useCount;
      if (--*useCount == 0)
      {
        if (WX_AR_IS_ALLOC(traits))
        {
          delete [] data;
          delete useCount;
        }
      }
      data = ad.data;
      useCount = ad.useCount;
      WX_AR_CLEAR_ALLOC(traits);
      return *this;
    }

    unsigned nelems;
    T *data;
    mutable int *useCount;
    unsigned traits;
};

#endif // __wxarraydata__
