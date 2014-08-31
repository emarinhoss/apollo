// warpx lib includes
#include <wxpetscindexer.h>

// std includes
#include <vector>

WxPetscIndexer::WxPetscIndexer(const WxSplitBox<int>& decomp)
  : idecomp(decomp) 
{
  unsigned nr = decomp.nboxes();
  offsets.push_back(0);
  WxIndexer<> tidx(decomp.getBox(0));
  idxrs.push_back(tidx);
  for (unsigned i=1; i<nr; ++i) 
  {
    WxIndexer<> ltidx(decomp.getBox(i));
    offsets.push_back( offsets[i-1] + decomp.getBox(i-1).area() );
    idxrs.push_back(ltidx);
  }
}
