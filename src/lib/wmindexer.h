#ifndef WmIndexer_h
#define WmIndexer_h

// std includes
#include <cassert>
#include <vector>

// WarpX includes
#include "wxrange.h"

/** \addtogroup lib 
 *  @{
 */

#define _WX_COL_MAJOR_ORDER 1
#define _WX_ROW_MAJOR_ORDER 2

// maxium rank array that can be indexed
const unsigned __max_idx_size=WxRange::max_dims;

///template<typename T> class WxArray;

/**
 * WmIndexer provides a way to map an element of an n-dimensional
 * index set into a single integer.
 */
template<int TYPE = _WX_COL_MAJOR_ORDER>
class WmIndexer
{
public:
	
	// this chummy-ness is presently required to make array slicing
	// work. Maybe fix this later.
	//template <typename T> friend class WxArray;
	
	/**
	 * The default constructor creates an empty WmIndexer object. This
	 * should be used only as a place holder.
	 */
	WmIndexer();
	
	/**
	 * Constructs indexer for a given range.
	 *
	 * @param r Range over which indexer is required
	 */
	WmIndexer(const WxRange& r);
	WmIndexer(unsigned rank, int *ai, const WxRange& r);
	~WmIndexer();
	
	WmIndexer(const WmIndexer& idx);
	WmIndexer& operator=(const WmIndexer& idx);
	
	/**
	 * Rank of indexer
	 */
	unsigned rank() const {
		return _rank;
	}
	
	/**
	 * Range of set indexed
	 */
	const WxRange& range() const {
		return _r;
	}
	
	/**
	 * Index 1D array
	 */
	int index(int k1) const;
	
	/**
	 * Index 2D array
	 */
	int index(int k1, int k2) const;
	
	/**
	 * Index 3D array
	 */
	int index(int k1, int k2, int k3) const;
	
	/**
	 * Index 4D array
	 */
	int index(int k1, int k2, int k3, int k4) const;
	
	/**
	 * Index arbitrary dimensional array
	 */
	int index(const int *k) const;
	
	/**
	 * Return index location given a linear offset
	 *
	 * User responsible to ensure return array large enough given range rank.
	 *
	 * @param loc linear offset
	 * @param indicesReturn indices as an array placed here.
	 * @return also returns address of returned indices
	 */
	int * invIndex(unsigned loc, int * indicesReturn ) const;
	
private:
	unsigned _rank;
	int _ai[__max_idx_size+1];
	WxRange _r;
};

template<int TYPE>
WmIndexer<TYPE>::WmIndexer()
{
}

template<int TYPE>
WmIndexer<TYPE>::WmIndexer(const WmIndexer<TYPE>& idx)
: _rank(idx._rank), _r(idx._r)
{
	for (unsigned i=0; i<=_rank; ++i) _ai[i] = idx._ai[i];
}

template<int TYPE>
WmIndexer<TYPE>::WmIndexer(unsigned rank, int *ai, const WxRange& r)
: _rank(rank), _r(r)
{
	for (unsigned i=0; i<=_rank; ++i) _ai[i] = ai[i];
}

template<int TYPE>
WmIndexer<TYPE>::~WmIndexer()
{
}

template<int TYPE>
WmIndexer<TYPE>&
WmIndexer<TYPE>::operator=(const WmIndexer<TYPE>& idx)
{
	if (this==&idx)
		return *this;
	
	_rank = idx._rank;
	for (int i=0; i<=_rank; ++i) _ai[i] = idx._ai[i];
	_r = idx._r;
	
	return *this;
}

//
// Col major order indexers
//
template<>
inline
int
WmIndexer<_WX_COL_MAJOR_ORDER>::index(int k1) const {
#ifdef _DO_RANGE_CHECK_
	assert(_rank==1);
	assert((k1>=_r.lower(0)) && (k1<=_r.upper(0)));
#endif 
	return _ai[0]+k1;
}

template<>
inline
int
WmIndexer<_WX_COL_MAJOR_ORDER>::index(int k1, int k2) const
{
#ifdef _DO_RANGE_CHECK_
	assert(_rank==2);
	assert((k1>=_r.lower(0)) && (k1<=_r.upper(0)));
	assert((k2>=_r.lower(1)) && (k2<=_r.upper(1)));
#endif 
	return _ai[0]+k1+_ai[2]*k2;
}

template<>
inline
int
WmIndexer<_WX_COL_MAJOR_ORDER>::index(int k1, int k2, int k3) const {
#ifdef _DO_RANGE_CHECK_
	assert(_rank==3);
	assert((k1>=_r.lower(0)) && (k1<=_r.upper(0)));
	assert((k2>=_r.lower(1)) && (k2<=_r.upper(1)));
	assert((k3>=_r.lower(2)) && (k3<=_r.upper(2)));
#endif 
	return _ai[0]+k1+_ai[2]*k2+_ai[3]*k3;         
}

template<>
inline
int
WmIndexer<_WX_COL_MAJOR_ORDER>::index(int k1, int k2, int k3, int k4) const {
#ifdef _DO_RANGE_CHECK_
	assert(_rank==4);
	assert((k1>=_r.lower(0)) && (k1<=_r.upper(0)));
	assert((k2>=_r.lower(1)) && (k2<=_r.upper(1)));
	assert((k3>=_r.lower(2)) && (k3<=_r.upper(2)));
	assert((k4>=_r.lower(3)) && (k4<=_r.upper(3)));
#endif  
	return _ai[0]+k1+_ai[2]*k2+_ai[3]*k3+_ai[4]*k4;        
}

template<>
inline
int *
WmIndexer<_WX_COL_MAJOR_ORDER>::invIndex(unsigned loc, int * indicesReturn) const
{
	int n = loc;
	div_t qr;
	for (int i=(int)_rank-1; i>=0; --i) {
		qr = div(n, _ai[i+1]);
		indicesReturn[i] = qr.quot + _r.lower(i);
		n = qr.rem;
	}
	return indicesReturn;
}

//
// Row major order indexers
//
template<>
inline
int
WmIndexer<_WX_ROW_MAJOR_ORDER>::index(int k1) const {
#ifdef _DO_RANGE_CHECK_
	assert(_rank==1);
	assert((k1>=_r.lower(0)) && (k1<=_r.upper(0)));
#endif 
	return _ai[0]+k1;
}

template<>
inline
int
WmIndexer<_WX_ROW_MAJOR_ORDER>::index(int k1, int k2) const
{
#ifdef _DO_RANGE_CHECK_
	assert(_rank==2);
	assert((k1>=_r.lower(0)) && (k1<=_r.upper(0)));
	assert((k2>=_r.lower(1)) && (k2<=_r.upper(1)));
#endif 
	return _ai[0]+_ai[1]*k1+k2;
}

template<>
inline
int
WmIndexer<_WX_ROW_MAJOR_ORDER>::index(int k1, int k2, int k3) const {
#ifdef _DO_RANGE_CHECK_
	assert(_rank==3);
	assert((k1>=_r.lower(0)) && (k1<=_r.upper(0)));
	assert((k2>=_r.lower(1)) && (k2<=_r.upper(1)));
	assert((k3>=_r.lower(2)) && (k3<=_r.upper(2)));
#endif 
	return _ai[0]+_ai[1]*k1+_ai[2]*k2+k3;         
}

template<>
inline
int
WmIndexer<_WX_ROW_MAJOR_ORDER>::index(int k1, int k2, int k3, int k4) const {
#ifdef _DO_RANGE_CHECK_
	assert(_rank==4);
	assert((k1>=_r.lower(0)) && (k1<=_r.upper(0)));
	assert((k2>=_r.lower(1)) && (k2<=_r.upper(1)));
	assert((k3>=_r.lower(2)) && (k3<=_r.upper(2)));
	assert((k4>=_r.lower(3)) && (k4<=_r.upper(3)));
#endif  
	return _ai[0]+_ai[1]*k1+_ai[2]*k2+_ai[3]*k3+k4;
}

template<>
inline
int *
WmIndexer<_WX_ROW_MAJOR_ORDER>::invIndex(unsigned loc, int * indicesReturn) const
{
	int n = loc;
	div_t qr;
	for (unsigned i=0; i<_rank; ++i) {
		qr = div(n, _ai[i+1]);
		indicesReturn[i] = qr.quot + _r.lower(i);
		n = qr.rem;
	}
	return indicesReturn;
}

template<int TYPE>
inline
int 
WmIndexer<TYPE>::index(const int *k) const {
#ifdef _DO_RANGE_CHECK_
	for (unsigned i=0; i<rank(); ++i)
		assert(k[i]>=_r.lower(i) && k[i]<=_r.upper(i));
#endif
	int sum=_ai[0];
	for (unsigned i=1; i<=_rank; ++i) sum += _ai[i]*k[i-1];
	return sum;
}

/** @}*/
#endif // WmIndexer_h
