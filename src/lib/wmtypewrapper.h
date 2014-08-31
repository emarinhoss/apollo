/*
 *  wmtypewrapper.h
 *  warpm Xcode project
 *
 *  provides utility to use std::type_info as key in associative containers
 *
 *  Created by Noah Reddell on 4/6/12.
 *  Copyright 2012 University of Washington. All rights reserved.
 *
 */
#ifndef wmtypewrapper_h
#define wmtypewrapper_h

#include <typeinfo>

/** \addtogroup lib 
 *  @{
 */


class WmTypeWrapper
{
public:
	WmTypeWrapper( const std::type_info &info ) : mInfo(info)
	{
	}
	
	// Requried to use as a key into a std::map.
	
	bool operator==( const WmTypeWrapper &other ) const
	{
		return mInfo.operator==( other.mInfo );
	}
	
	bool operator<( const WmTypeWrapper &other ) const
	{
		return mInfo.before( other.mInfo );
	}
	
	bool operator==( const std::type_info & other ) const
	{
		return mInfo.operator==( other );
	}
	
	
private:
	const std::type_info &mInfo;

};



#endif // wmtypewrapper_h