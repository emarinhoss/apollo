#ifndef wxsplitrange_h
#define wxsplitrange_h

// WarpX lib includes
#include <wxrange.h>

// std includes
#include <map>
#include <vector>
#include <list>
#include <cmath>

/**
 * WxSplitRange splits a range into smaller ranges in a cartesian manner.  Each sub-range face has a 1:1 mapping with neighbor sub-range faces.
 */
template <typename ID_type>
class WxSplitRange
{
public:
	
	typedef int empty_range_t;
	typedef std::vector< ID_type > SubRangesAssigned_C_ordering_t;
	
	/**
	 * Split supplied range into smaller number of ranges.
	 *
	 * Each sub-range is labelled using elements from the supplied list of ids. 
	 * Each id is used no more than once. All ids will not necessarily be used.
	 *
	 * @param range range to split
	 * @param id_list list of ids to give to sub-ranges. Ordering of assignment will be as done with C-style row-major counting. (incrementing highest dimension first)
	 */
	WxSplitRange( const WxRange& range, const std::list<ID_type> & id_list )  : _ndims(range.ndims()), _fullRange(range), _lowerPoints(range.ndims()), _upperPoints(range.ndims())
	{
		_numsubranges = _divideRange( id_list.size(), range, _lowerPoints, _upperPoints );
		// create map to subranges
		std::vector<int> subrangeIndex( _ndims, 0 );
		typename std::list<ID_type>::const_iterator id_itr = id_list.begin();
		do {
			_subrangeMap.insert( SubRangePair_t( *id_itr, subrangeIndex ) );
			_assignedSubRange_ids.push_back(*id_itr);
			++id_itr;
		} while (_advanceIndex( subrangeIndex ));		
	};
	
	/**
	 * Destory a split range
	 */
	virtual ~WxSplitRange() {};
	
	/**
	 * Dimension of range split
	 */
	unsigned ndims() const {
		return _fullRange.ndims();
	}
	
	/**
	 * Number of sub-ranges in split
	 */
	unsigned numRanges() const {
		return _numsubranges;
	}
	
	
	/**
	 * Return the global range that was split.
	 */
	const WxRange & getGlobalRange() const {
		return _fullRange;
	}
	
	/**
	 * Returns the assigned subrange ids for all of the assigned subranges
	 *
	 * Ordering of vector returned is as per standard C-style ordering or highest dimension subrange index varying the fastest.
	 */
	const SubRangesAssigned_C_ordering_t getAssignedSubRangeIDs() const
	{
		return _assignedSubRange_ids;
	}
	
	/**
	 * Return neighbor for a given sub-range and side
	 *
	 * @param id sub-rangefor which neighbor information is needed
	 * @param dim dimension in which the adjacent side face lies
	 * @param lower_or_upper in which direction in dimension dim to look for neighbor
	 * @return neighboring sub-range id, or ??? if no neighbor
	 */
	/*  To be implemented later when we need it
	 ID_type getNeigh() {
      return neighHolder->getNeigh(n, *this, lp, up);
	 }
	 */
	
	/**
	 * Return a sub-range in the split.
	 *
	 * @param id sub-range number. Should be in [0..nranges)
	 * @return sub-range associated with 'id'
	 */
	WxRange getRange( const ID_type & id ) const
	{
		typename SubRangeMap_t::const_iterator itr = _subrangeMap.find(id);
		if (itr != _subrangeMap.end() )
		{
			Distribution_Position_t index = (*itr).second;
			std::vector<int> lower, upper;
			for (int i=0; i<_ndims; ++i) {
				lower.push_back( _lowerPoints[i][index[i]] );
				upper.push_back( _upperPoints[i][index[i]] );
			}
			return WxRange( _ndims, &lower[0], &upper[0] );
		}

		// not found: throw an exception
		empty_range_t e;
		throw e;
	}
	
	const std::vector<int> & getLowerPoints( int dim ) const
	{
		return _lowerPoints[dim];
	}

	const std::vector<int> & getUpperPoints( int dim ) const
	{
		return _upperPoints[dim];
	}
	
	/**
	 * modifies the split range such that in the dimension specified
	 * all the range coordinates are stretched in the positive direction
	 * by the supplied factor. (lower range point remains fixed.)
	 **/
	void stretchDimension( int dim, unsigned int stretchFactor)
	{
		std::vector<int> lowerExtension(_ndims, 0);
		std::vector<int> upperExtension(_ndims, 0);
		upperExtension[dim] = _fullRange.length(dim) * (stretchFactor-1);
		
		_fullRange = _fullRange.extend(&lowerExtension[0], &upperExtension[0]);
		
		int numPoints = _lowerPoints[dim].size();
		std::vector<int> newLower(numPoints);
		int baseline = _lowerPoints[dim][0]; //starting point of lowest range remains the same
		
		for (int i=0; i<numPoints; ++i) {
			newLower[i] = (_lowerPoints[dim][i] - baseline )* stretchFactor + baseline;
			_upperPoints[dim][i] = (_upperPoints[dim][i] - _lowerPoints[dim][i] + 1)*stretchFactor + newLower[i] - 1;
		}
		_lowerPoints[dim] = newLower;
	}
	
	/**
	 * modifies the split range such that the number of dimensions is
	 * increased by one, covering the unplit range from lowerPoint to
	 * upperPoint.  Previous splitting in existing dimensions remains the same.
	 **/
	void addDimension( int lowerPoint, int upperPoint )
	{
		std::vector<int> newlowervec( 1, lowerPoint);
		std::vector<int> newuppervec( 1, upperPoint);

		_lowerPoints.push_back( newlowervec );
		_upperPoints.push_back( newuppervec );
		_ndims += 1;
		_fullRange = _fullRange.extDim(lowerPoint, upperPoint);

		// lastly have to add range index of 0 to the new subrange indices
		for( typename SubRangeMap_t::iterator itr = _subrangeMap.begin(); itr != _subrangeMap.end(); ++itr )
			(*itr).second.push_back(0);
	}
	
	
	/**
	 * Given a point, return which range it resides in
	 *
	 * @param coords coordinates of the point in question
	 * @return id of the sub-range in which the point resides
	 */
	/* const ID_type & getRank( const std::vector<int> coords ) const
	{
		unsigned
		WxSplitRange<ID_type>::getRank(TYPE coords[]) const 
		{
			typename SubRangeMap_t::const_iterator i, iend;
			iend = _data->rangeMap.end();
			for (i = _data->rangeMap.begin(); i != iend; ++i) {
				if ((*i).second.contains(coords)) {
					return (*i).first;
				}
			}
			//  Doesn't exist in the decomp: throw a tde.
			WxExcept wxe("WxSplitRange::getRank: ");
			wxe << "Point " << coords[0];
			for (unsigned i = 1; i < _ndims; ++i)
				wxe << "," << coords[i];
			wxe << " is not in decomposition.";
			throw wxe;
		}	
	} */
	
private:
	
	// types of ID_type -> WxRange map and pairs
	typedef std::vector<int> Distribution_Position_t;  // split coordinates
	typedef std::map<ID_type, Distribution_Position_t > SubRangeMap_t;
	typedef std::pair<ID_type, Distribution_Position_t > SubRangePair_t;
	
	typedef std::vector< std::vector<int> > Distribution_Storage_t;

	
	unsigned _ndims;
	WxRange _fullRange;
	Distribution_Storage_t _lowerPoints;
	Distribution_Storage_t _upperPoints;
	SubRangeMap_t _subrangeMap;
	int _numsubranges;
	SubRangesAssigned_C_ordering_t _assignedSubRange_ids;
	
	// TODO: certainly more can be done here, but this is generally working now. In the future, more can be done to maximize subrange volume to surface area ratio and maximally use subranges available
	int _divideRange( unsigned maxSubRangesDesired, const WxRange& decompRange, Distribution_Storage_t & lowerPoints, Distribution_Storage_t & upperPoints )
	{
		// first-pass attempt at dividing the range up into NxM sub-ranges.  arrangement of subranges structured to be compatible with Global Arrays' irreg config call
		
		int ndims = decompRange.ndims();
		//float averageVolume = decompRange.size() / (float) maxSubRangesDesired;
		float uniformNumSplits_f = std::pow( maxSubRangesDesired, 1.0 / ndims );
		int numSplitsInt = std::floor( uniformNumSplits_f );  //nominal number of sub-ranges (splits+1) in each dimension
		
		int numSubRangesActuallyAssigned = 1;
		
		for (int i = 0; i < ndims; ++i) {
			int numThisDim = 0;
			int nominalLength = (float) decompRange.length(i) / numSplitsInt;
			if ( (i == ndims-1) ) // if we can afford some more splits in this last dim, do so
			{
				int remainingSplitsThisDim = maxSubRangesDesired / numSubRangesActuallyAssigned;
				nominalLength = (float) decompRange.length(i) / ((float) remainingSplitsThisDim );
			}
			int lowerPos = decompRange.lower(i);
			do {
				lowerPoints[i].push_back( lowerPos );
				int upperPos = std::min( lowerPos + nominalLength -1, decompRange.upper(i) );
				if (i == ndims-1)
					if ( (numThisDim+2)*numSubRangesActuallyAssigned > maxSubRangesDesired ) // special case to prevent exceeding maxSubRangesDesired
						upperPos = decompRange.upper(i);// this should be the last split added, so extend upper pos to end
						

				upperPoints[i].push_back( upperPos );
				lowerPos = upperPos+1;
				++numThisDim;
			} while (lowerPos <= decompRange.upper(i) );
			numSubRangesActuallyAssigned *= numThisDim;
		}
		return numSubRangesActuallyAssigned;
	};
	
	bool _advanceIndex( Distribution_Position_t & subrangeIndex, int dim = 0 ) const
	{
		bool higherResult;
		if (dim < _ndims -1 )
		{	higherResult = _advanceIndex( subrangeIndex, dim+1 ); // recursively advance
			if (higherResult)
				return true;
		}
		// reaching here, we need to advance this dim or return false if we can't
			
		int nextIndex = subrangeIndex[dim] + 1;
		if (nextIndex == _lowerPoints[dim].size() )
			return false; // nowhere further to go
		
		for (int i=dim+1; i < _ndims; ++i)
			subrangeIndex[i] = 0; // reset all higher indices to zero
		
		subrangeIndex[dim] = nextIndex; // advance this index
		return true;
	};
};

#endif //  wxsplitrange_h
