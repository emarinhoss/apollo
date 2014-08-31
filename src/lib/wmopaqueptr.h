#ifndef wmopaqueptr_h
#define wmopaqueptr_h
//
//  wmopaqueptr.h
//  warpm Xcode project
//
//  Created by Noah Reddell on 8/9/11.
//  Copyright 2011 University of Washington. All rights reserved.
//

#include <wmtypewrapper.h>

/** \addtogroup lib 
 *  @{
 */

/**
 * Opaque Pointer object that carries the referenced data type, but does not expose
 * this payload type as part of the poiter type.
 **/
class WmOpaquePtr
{
public:
	WmOpaquePtr( void * data, const std::type_info & type): _data(data), _type(type) {};
	
	
	const WmTypeWrapper getType() const {return _type;};
	void * getDataPtr() const { return _data; };
	
private:
	WmOpaquePtr(); //no reason to allow an undefined instance, so block default constructor
	
	void * _data;
	const WmTypeWrapper _type;
};

/**
 * Opaque Pointer object that carries the referenced data type, but does not expose
 * this payload type as part of the poiter type.
 * 
 * Points to constant payload.
 **/
class WmConstOpaquePtr
{
public:
	WmConstOpaquePtr( const void * data, const std::type_info & type): _data(data), _type(type) {};
	
	
	const WmTypeWrapper getType() const {return _type;};
	const void * getDataPtr() const { return _data; };
	
private:
	WmConstOpaquePtr(); //no reason to allow an undefined instance, so block default constructor
	
	const void * _data;
	const WmTypeWrapper _type;
};


/** @}*/
#endif //__wmopaqueptr_h
