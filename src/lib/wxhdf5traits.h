#ifndef __wxhdf5traits__
#define __wxhdf5traits__

// hdf5
#define H5_USE_16_API
#include <hdf5.h>

/**
 * Traits class for HDF5
 */
template <typename T> struct WxHdf5Traits;

/**
 * HDF5 traits class for char
 */
template <> struct WxHdf5Traits<char> {
    static hid_t hdf5Type() {
        return H5T_NATIVE_CHAR;
    }
};

/**
 * HDF5 traits class for unsigned char
 */
template <> struct WxHdf5Traits<unsigned char> {
    static hid_t hdf5Type() {
        return H5T_NATIVE_UCHAR;
    }
};

/**
 * HDF5 traits class for short
 */
template <> struct WxHdf5Traits<short> {
    static hid_t hdf5Type() {
        return H5T_NATIVE_SHORT;
    }
};

/**
 * HDF5 traits class for unsigned short
 */
template <> struct WxHdf5Traits<unsigned short> {
    static hid_t hdf5Type() {
        return H5T_NATIVE_USHORT;
    }
};

template <> struct WxHdf5Traits<int> {
    static hid_t hdf5Type() {
        return H5T_NATIVE_INT;
    }
};

template <> struct WxHdf5Traits<unsigned int> {
    static hid_t hdf5Type() {
        return H5T_NATIVE_UINT;
    }
};

template <> struct WxHdf5Traits<long> {
    static hid_t hdf5Type() {
        return H5T_NATIVE_LONG;
    }
};

template <> struct WxHdf5Traits<unsigned long> {
    static hid_t hdf5Type() {
        return H5T_NATIVE_ULONG;
    }
};

template <> struct WxHdf5Traits<float> {
    static hid_t hdf5Type() {
        return H5T_NATIVE_FLOAT;
    }
};

template <> struct WxHdf5Traits<double> {
    static hid_t hdf5Type() {
        return H5T_NATIVE_DOUBLE;
    }
};

template <> struct WxHdf5Traits<long double>  {
    static hid_t hdf5Type() {
        return H5T_NATIVE_LDOUBLE;
    }
};

template <> struct WxHdf5Traits<long long> {
    static hid_t hdf5Type() {
        return H5T_NATIVE_LLONG;
    }
};

#endif // __wxhdf5traits__
