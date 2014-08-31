#include "wxhdf5traits.h"

#include <map>
#include <boost/assign.hpp>


/* here we populate the static mapping of native C++ data types to the correspoiding HDF5
 * data type identifier.
 */

WxHdf5Traits::hdf5TypeMap_t WxHdf5Traits::hdf5TypeMap = boost::assign::map_list_of 
( WmTypeWrapper( typeid(char) ), H5T_NATIVE_CHAR )
( WmTypeWrapper( typeid(unsigned char) ), H5T_NATIVE_UCHAR )
( WmTypeWrapper( typeid(short) ), H5T_NATIVE_SHORT )
( WmTypeWrapper( typeid(unsigned short) ), H5T_NATIVE_USHORT )
( WmTypeWrapper( typeid(int) ), H5T_NATIVE_INT )
( WmTypeWrapper( typeid(unsigned int) ), H5T_NATIVE_UINT )
( WmTypeWrapper( typeid(long) ), H5T_NATIVE_LONG )
( WmTypeWrapper( typeid(unsigned long) ), H5T_NATIVE_ULONG )
( WmTypeWrapper( typeid(float) ), H5T_NATIVE_FLOAT )
( WmTypeWrapper( typeid(double) ), H5T_NATIVE_DOUBLE )
( WmTypeWrapper( typeid(long double) ), H5T_NATIVE_LDOUBLE )
( WmTypeWrapper( typeid(long long) ), H5T_NATIVE_LLONG );
