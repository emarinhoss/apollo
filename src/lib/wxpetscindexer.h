#ifndef __wxpetscindexer__h__
#define __wxpetscindexer__h__

// warpx lib includes
#include <wxsplitbox.h>
#include <wxindexer.h>

//  std includes
#include <vector>

/**
 * Indexing routines designed to align data closely with the PETSc
 * data system
 */
class WxPetscIndexer {
  public:
/**
 * Create new indexer based on given decomposition
 *
 * @param decomp the decomposition
 */
    WxPetscIndexer(const WxSplitBox<int>& decomp);

/**
 * Return linear index given a two dimensional (i,j) index
 *
 * @param i the i'th index
 * @param j the j'th index
 * @return Linear index into array
 */
    int index(int i, int j) {
      int idx[2];
      idx[0] = i;
      idx[1] = j;
      unsigned rank = idecomp.getRank(idx);
      return offsets[rank] + idxrs[rank].index(idx);
    }

/**
 * Return (i,j) given a linear index
 *
 * @param n linear index
 * @return (i,j) index corresponding to n
 */
    const int* invIndex(int n) {
      std::vector<int>::iterator oitr = offsets.end();
      std::vector<WxIndexer<> >::iterator iitr = idxrs.end();
      do {
        --oitr;
        --iitr;
      } while ( n < *oitr );
      
      return (*iitr).invIndex(n - *oitr);
    }

  private:
/** Array of offsets */
    std::vector<int> offsets;
/** Array of box indexers */
    std::vector<WxIndexer<> > idxrs;
/** Decomp on which this indexer is based */
    WxSplitBox<int> idecomp;
/** Coordinates this indexer is currently representing */
    int indices[2];
};

#endif // __wxpetscindexer__h__

