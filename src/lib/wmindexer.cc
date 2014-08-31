#include "wmindexer.h"

template<>
WmIndexer<_WX_COL_MAJOR_ORDER>::WmIndexer(const WxRange &r)
: _rank(r.ndims()), _r(r)
{
	// set a_1, ... a_N
	_ai[1] = 1;
	for (unsigned i=2; i<=_rank; ++i) 
		_ai[i] = _ai[i-1]*r.length(i-2);
	
	// set a_0
	int sum = 0;
	for (unsigned i=1; i<=_rank; ++i) sum += _ai[i]*r.lower(i-1);
	_ai[0] = -sum;
}

template<>
WmIndexer<_WX_ROW_MAJOR_ORDER>::WmIndexer(const WxRange &r)
: _rank(r.ndims()), _r(r)
{
	// set a_1, ... a_N
	_ai[_rank] = 1;
	for (int i=_rank-1; i>=1; --i) 
		_ai[i] = _ai[i+1]*r.length(i);
	
	// set a_0
	int sum = 0;
	for (int i=1; i<=_rank; ++i) sum += _ai[i]*r.lower(i-1);
	_ai[0] = -sum;
}

// instantiations
template class WmIndexer<_WX_COL_MAJOR_ORDER>;
template class WmIndexer<_WX_ROW_MAJOR_ORDER>;
