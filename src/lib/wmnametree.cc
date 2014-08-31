/*
 *  wmnametree.cc
 *  warpm Xcode project
 *
 *  Created by Noah Reddell on 3/30/11.
 *  Copyright 2011 University of Washington. All rights reserved.
 *
 */

#include "wmnametree.h"

template <typename V>
WmNameTree<V> WmNameTree<V>::find_children( const std::string & baseNodeName ) const
{
	using namespace std;
	WmNameTree<V> childMatches;
	
	string prefix( baseNodeName+"." ); // parent prefix plus the '.' delimiter
	
	typename map<string, V >::const_iterator lastMatch = childMatches.begin(); 
	typename map<string, V >::const_iterator endpoint = this->end(); 
	// later this could be implemented more efficiently by using a different sort function in the parent map class
	// for now, we will iterate through each member, checking for components that match the base node prefix
	for (typename map<string, V >::const_iterator it = this->begin(); it != endpoint; ++it)
	{
		if (! it->first.compare( 0, prefix.size(), prefix ) )
		{
			//both strings begin with the same prefix, so add this element to the list of matches - without prefix
			//we use lastMatch as a hint for insertion.  Which should be dead-on, because we are just copying elemens over in order
			childMatches.insert( lastMatch++, pair<string, V>( it.first.substr(prefix.size()) /*remove prefix*/ , it.second )  );
			// TODO: this lastMatch hint is likely no longer good now that I am stripping prefix.  It all comes down to using the defaulet
			// sorting function for strings.  I should probably implement my own for performance.
		}
	}
	
	
	return childMatches;
	
}


template <typename V>
unsigned WmNameTree<V>::num_children( const std::string & baseNodeName ) const
{
	//this will also include grandchildren, etc.
	WmNameTree<V> childMatches = find_children(baseNodeName);
	return childMatches.size();
	
}