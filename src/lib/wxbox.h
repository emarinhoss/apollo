#ifndef __wxbox__h__
#define __wxbox__h__

// WarpX includes
#include "wxobject.h"
#include "wxcryptset.h"

// std includes
#include <cassert>
#include <iostream>
#include <vector>

/**
 * WxBox represents a n-dimensional box specified by lower bounds and
 * upper bounds.
 */
template<typename TYPE>
class WxBox : public WxObject
{
  public:

    WxBox() 
      : _ndims(0) {
    }

/**
 * Create a 'ndims' dimensional box
 */
    WxBox(unsigned ndims);

/**
 * Constucts a box with given 'lower' and 'upper' bounds.
 *
 * @param ndims Dimensionality of the box
 * @param lower Array of lower corner coordinates
 * @param upper Array of upper corner coordinates
 */
    WxBox(unsigned ndims, TYPE *lower, TYPE *upper);

/**
 * Constucts a box given lengths of each side.
 *
 * @param ndims Dimensionality of the box
 * @param length Array of upper corner coordinates
 */
    WxBox(unsigned ndims, TYPE *length);

/**
 * Rank-1 ctor
 */
    WxBox(TYPE s1, TYPE e1);
/**
 * Rank-2 ctor
 */
    WxBox(TYPE s1, TYPE e1, TYPE s2, TYPE e2);

/**
 * Rank-3 ctor
 */
    WxBox(TYPE s1, TYPE e1, TYPE s2, TYPE e2, TYPE s3, TYPE e3);

/**
 * Rank-4 ctor
 */
    WxBox(TYPE s1, TYPE e1, TYPE s2, TYPE e2, TYPE s3, TYPE e3, TYPE s4, TYPE e4);

/**
 * Copy constructor and assignment operators
 */
    WxBox(const WxBox& b);
    WxBox& operator=(const WxBox& b);

    virtual ~WxBox();

/**
 * Constructs a box specified in a cryptset. The cryptset should have
 * variables
 *   
 * ndims = int
 * lower = [TYPE, ...]
 * upper = [TYPE, ...]
 *
 * @param wxc WxCryptSet to initialize box with.
 */
    virtual void setup(const WxCryptSet& wxc);

/**
 * Dimensionality of box
 *
 * @return dimensionality of box
 */
    unsigned ndims() const {
      return _ndims;
    }

/**
 * Lower bound along dimension 'dim'.
 *
 * @param dim Direction
 * @return lower bound along dimension dim
 */
    TYPE lower(unsigned dim) const { 
      return _lower[dim]; 
    }

/**
 * Set lower bound along dimension 'dim'.
 *
 * @param dim Direction
 * @return lower bound along dimension dim
 */
    void lower(unsigned dim, TYPE low) { 
      _lower[dim] = low;
      _length[dim] = _upper[dim] - _lower[dim];
    }

/**
 * Upper bound along dimension 'dim'.
 *
 * @param dim Direction
 * @return Upper bound along dimension dim
 */
    TYPE upper(unsigned dim) const { 
      return _upper[dim]; 
    }

/**
 * Set upper bound along dimension 'dim'.
 *
 * @param dim Direction
 * @return Upper bound along dimension dim
 */
    void upper(unsigned dim, TYPE upp) { 
      _upper[dim] = upp;
      _length[dim] = _upper[dim] - _lower[dim];
    }

/**
 * Determine if box is empty
 *
 * @return true of box is empty, false otherwise
 */
    bool isEmpty() const;
    
/**
 * Length of edge along dimension 'dim'.
 *
 * @param dim Direction
 * @return Lenght of side along dimension dim
 */
    TYPE length(unsigned dim) const {
      return _length[dim];
    }

/**
 * Area of box
 *
 * @return area of box
 */
    TYPE area() const;

/**
 * Area of box
 *
 * @return area of box
 */
    TYPE size() const {
      return area();
    }

/**
 * Computes intersection of this box with supplied box
 *
 * @param box to intersect with
 * @return intersection box
 */
    WxBox<TYPE> intersect(const WxBox<TYPE>& box) const;


/**
* Computes intersection of this box with supplied box
*
* @param box to intersect with
* @return intersection box.
* It only checks the transverse dimensions.  Perpendicular
* dimension is set to the this-> box boundary condition.
*/

    WxBox<TYPE> intersectmultiblock(const WxBox<TYPE>& box, unsigned dir, int edge1, int edge2, bool recv, int npad) const;

/**
 * Returns a box which is extended by the given amount along each side
 * in each dimension
 *
 * @param low extend region along lower edge in direction dim is low[dim]
 * @param upp extend region along upper edge in direction dim is upp[dim]
 * @return extended box
 */
    WxBox<TYPE> extend(const TYPE low[], const TYPE upp[]) const;

/**
 * Returns a new box which has one greater dimension that this
 * one.
 *
 * @param low lower coordinate for new dimension
 * @param upp upper coordinate for new dimension
 * @return box of one higher dimension. This box has the new dimension as its
 *         first one.
 */
    WxBox<TYPE> extDim(TYPE low, TYPE upp) const;

/**
 * Check if this box is equal to one supplied
 *
 * @param b Box to compare with
 * @return true if the boxes are same, or false otherwise
 */
    bool operator==(const WxBox<TYPE>& b) const;

/**
 * Check if the box contains the given point
 *
 * @param coord Point to check
 * @return true if box contains this point, false otherwise
 */
    bool contains(TYPE coord[]) const;

  private:
    unsigned _ndims;
    TYPE _lower[16], _upper[16], _length[16];
};

/**
 * Returns true if boxes are same
 */
template<typename TYPE>
bool operator==(const WxBox<TYPE> &wa, const WxBox<TYPE>& wb);

#endif // __wxbox__h__
