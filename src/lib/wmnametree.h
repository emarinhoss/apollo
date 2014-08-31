/*
 *  wmnametree.h
 *  warpm Xcode project
 *
 *  Created by Noah Reddell on 3/30/11.
 *  Copyright 2011 University of Washington. All rights reserved.
 *
 */

#ifndef nametree_hh_
#define nametree_hh_

#include <map>
#include <string>

/** \addtogroup lib 
 *  @{
 */

/** Provides hierarchical container for T values according to unique string keys. 
 *
 * Provides map class functionality where the unique keys are strings which are expected to be hierarchical.
 * Hierarchy is determined by '.' delemiters in the string key.
 * Ex: density, velocity.x, velocity.y, current.ion.x, current.electron.y are all valid keys.
 * From above, velocity.x and velocity.y are both children of the velocity node.  The parent node may or may not
 * be an actual member of the container.  Above, note no element 'velocity', just children of velocity.
 **/
template<typename V>
class WmNameTree : public std::map<std::string, V >  {
	
public:
	
	/*! Basic constructor creates empty container for elements of type V with unique string keys.
	 *
	 */
	WmNameTree<V>() : std::map<std::string, V >() {};
	
	/*! Return all children (and grandchildren, etc.) of a baseNode.  baseNodeName prefix is stripped.
	 *
	 * @param baseNodeName The string identifying the base node
	 * @return A new tree holding only the children.
	 */
	WmNameTree<V> find_children( const std::string & baseNodeName ) const;
	
	/*! Return count of all children (and grandchildren, etc.) of a baseNode.
	 *
	 * @param baseNodeName The string identifying the base node
	 * @return Number of children and grandchildren, etc.
	 */
	unsigned num_children( const std::string & baseNodeName ) const;
	
	
	/* *** No Member variables!!!!! ***  Deriving from a std::map is not recommended, because its destructor is not virtual
	 I chose to do it anyway, because I do not see a need for new members variables while the std::map member functions
	 still make sense */
	
	
private:
	
};



/** @}*/
#endif //nametree_hh_