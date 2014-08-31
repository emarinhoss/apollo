#ifndef __wxgridbox__h__
#define __wxgridbox__h__

// WarpX includes
#include "wxrange.h"
#include "wxbox.h"
#include "wxcryptset.h"
#include "wxgridrange.h"

// std includes
#include <vector>

/**
 * WxGridBox extends WxBox class so that it can be divided into
 * cells. Thus it represents a gridded rectangular domain of arbitrary
 * dimension.
 */
template<typename REAL>
class WxGridBox : public WxBox<REAL>, public WxGridRange
{
  public:

/**
 * Does nothing: only a place holder
 */
    WxGridBox() {
    }
    
/**
 * Create a 'ndims' dimensional grid box
 */
    WxGridBox(unsigned ndims);

/**
 * Constructs a new gridbox.
 *
 * @param dims Dimension of box
 * @param lower Coordinates of lower corner of box
 * @param upper Coordinates of upper corner of box
 * @param r Range object for box interior
 * @param npad Ghost cells on each side of box
 */
    WxGridBox(unsigned dims, REAL *lower, REAL *upper, const WxRange& r, unsigned npad=0);

    WxGridBox(const WxGridBox& gb);
    WxGridBox& operator=(const WxGridBox& gb);

    virtual ~WxGridBox() {}

/**
 * Constructs grid-box from supplied cryptset. The cryptset should
 * have variables
 *
 * ndims = int
 * lower = [REAL, ...]
 * upper = [REAL, ...]
 * npad = int
 * cells = [int, ...]
 *
 * @param wxc WxCryptSet to initialize gridbox with.
 */
    virtual void setup(const WxCryptSet& wxc);

/**
 * Return manual decomp associated with the grid
 *
 * @return manual decomp. Empty list is returned if there was no
 * decomp specified.
 */
    std::vector<unsigned> getManualDecompSpec() {
      return _decomp;
    }

/**
 * Dimension of grid
 */
    unsigned ndims() const {
      return WxBox<REAL>::ndims();
    }    

/**
 * Cell size along dimension 'dim'
 */
    REAL dx(unsigned dim) const {
      return this->length(dim)/(REAL) this->range().length(dim);
    }

    REAL femdx(unsigned dim) const {
      return this->length(dim)/((REAL) this->range().length(dim)-1);
    }

/**
 * Coordinates of cell center
 *
 * @param idx Index of cell
 * @param x On output the coordinate of the cell-center
 */
    void coord(const int *idx, REAL *x) const;

/**
 * Coordinates of cell center for 1D block
 *
 * @param i Index of cell
 * @param x On output the coordinate of the cell-center
 */
    void coord(int i, REAL *x) const;

/**
 * Coordinates of cell center for 2D block
 *
 * @param i,j Index of cell
 * @param x On output the coordinate of the cell-center
 */
    void coord(int i, int j, REAL *x) const;

/**
 * Coordinates of cell center for 3D block
 *
 * @param i,j,k Index of cell
 * @param x On output the coordinate of the cell-center
 */
    void coord(int i, int j, int k, REAL *x) const;

/**
 * Coordinates of node
 *
 * @param idx Index of cell
 * @param x On output the coordinate of the node
 */
    void nodeCoord(const int *idx, REAL *x) const;
    void femNodeCoord(const int *idx, REAL *x) const;

/**
 * Coordinates of node for 1D block
 *
 * @param i Index of cell
 * @param x On output the coordinate of the node
 */
    void nodeCoord(int i, REAL *x) const;
    void femNodeCoord(int i, REAL *x) const;

/**
 * Coordinates of node for 2D block
 *
 * @param i,j Index of cell
 * @param x On output the coordinate of the node
 */
    void nodeCoord(int i, int j, REAL *x) const;

/**
 * Coordinates of node for 3D block
 *
 * @param i,j,k Index of cell
 * @param x On output the coordinate of the node
 */
    void nodeCoord(int i, int j, int k, REAL *x) const;

/**
 * Return true if domain in direction 'dir' is periodic
 *
 * @param dir Direction to check
 * @return true if direction is periodic, false otherwise
 */
    bool isDirPeriodic(unsigned dir) const {
      return _isDirPeriodic[dir];
    }

/**
 * Return array indicating periodic directions
 *
 * @return array indicating periodic directions
 */
    bool* periodicDirs() {
      return _isDirPeriodic;
    }

  private:
// indicates which directions are periodic
    bool _isDirPeriodic[16];
// manual decomp for this grid
    std::vector<unsigned> _decomp;
};

template<typename REAL>
inline
void 
WxGridBox<REAL>::coord(const int *idx, REAL *x) const 
{
  for (unsigned i=0; i<this->ndims(); ++i)
    x[i] = this->lower(i) +
      this->dx(i)*(idx[i]-this->lowerIdx(i) + 0.5);
}

template<typename REAL>
inline
void 
WxGridBox<REAL>::coord(int i, REAL *x) const 
{
  x[0] = this->lower(0) +
    this->dx(0)*(i - this->lowerIdx(0) + 0.5);
}

template<typename REAL>
inline
void 
WxGridBox<REAL>::coord(int i, int j, REAL *x) const 
{
  x[0] = this->lower(0) +
    this->dx(0)*(i - this->lowerIdx(0) + 0.5);
  x[1] = this->lower(1) +
    this->dx(1)*(j - this->lowerIdx(1) + 0.5);
}

template<typename REAL>
inline
void 
WxGridBox<REAL>::coord(int i, int j, int k, REAL *x) const 
{
  x[0] = this->lower(0) +
    this->dx(0)*(i - this->lowerIdx(0) + 0.5);
  x[1] = this->lower(1) +
    this->dx(1)*(j - this->lowerIdx(1) + 0.5);
  x[2] = this->lower(2) +
    this->dx(2)*(k - this->lowerIdx(2) + 0.5);
}


template<typename REAL>
inline
void
WxGridBox<REAL>::nodeCoord(const int *idx, REAL *x) const 
{
  for (unsigned i=0; i<this->ndims(); ++i)
    x[i] = this->lower(i) +
      this->dx(i)*(idx[i]-this->lowerIdx(i));
}

template<typename REAL>
inline
void 
WxGridBox<REAL>::nodeCoord(int i, REAL *x) const 
{
  x[0] = this->lower(0) +
    this->dx(0)*(i - this->lowerIdx(0));
}

template<typename REAL>
inline
void
WxGridBox<REAL>::femNodeCoord(const int *idx, REAL *x) const
{
  for (unsigned i=0; i<this->ndims(); ++i)
    x[i] = this->lower(i) +
      this->femdx(i)*(idx[i]-this->lowerIdx(i));
}

template<typename REAL>
inline
void
WxGridBox<REAL>::femNodeCoord(int i, REAL *x) const
{
  x[0] = this->lower(0) +
    this->femdx(0)*(i - this->lowerIdx(0));
}

template<typename REAL>
inline
void 
WxGridBox<REAL>::nodeCoord(int i, int j, REAL *x) const 
{
  x[0] = this->lower(0) +
    this->dx(0)*(i - this->lowerIdx(0));
  x[1] = this->lower(1) +
    this->dx(1)*(j - this->lowerIdx(1));
}

template<typename REAL>
inline
void
WxGridBox<REAL>::nodeCoord(int i, int j, int k, REAL *x) const 
{
  x[0] = this->lower(0) +
    this->dx(0)*(i - this->lowerIdx(0));
  x[1] = this->lower(1) +
    this->dx(1)*(j - this->lowerIdx(1));
  x[2] = this->lower(2) +
    this->dx(2)*(k - this->lowerIdx(2));
}

#endif // __wxgridbox__h__
