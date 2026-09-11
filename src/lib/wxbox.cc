#include "wxbox.h"

template<typename TYPE>
WxBox<TYPE>::WxBox(unsigned ndims)
  : _ndims(ndims)
{
  for (unsigned i=0; i<_ndims; ++i) 
  {
    _lower[i] = 0;
    _upper[i] = 0;
    _length[i] = 0;
  }    
}

template<typename TYPE>
WxBox<TYPE>::WxBox(unsigned ndims, TYPE *lower, TYPE *upper)
  : _ndims(ndims)
{
  for (unsigned i=0; i<_ndims; ++i) 
  {
    _lower[i] = lower[i];
    _upper[i] = upper[i];
    _length[i] = _upper[i]-_lower[i];
  }
}

template<typename TYPE>
WxBox<TYPE>::WxBox(unsigned ndims, TYPE *length)
  : _ndims(ndims)
{
  for (unsigned i=0; i<_ndims; ++i) 
  {
    _lower[i] = 0;
    _upper[i] = length[i];
    _length[i] = _upper[i]-_lower[i];
  }
}

template<typename TYPE>
WxBox<TYPE>::WxBox(TYPE s1, TYPE e1)
  : _ndims(1)
{
  _lower[0] = s1; _upper[0] = e1;
  for (unsigned i=0; i<_ndims; ++i) 
    _length[i] = _upper[i]-_lower[i];
}

template<typename TYPE>
WxBox<TYPE>::WxBox(TYPE s1, TYPE e1, TYPE s2, TYPE e2)
  : _ndims(2)
{
  _lower[0] = s1; _upper[0] = e1;
  _lower[1] = s2; _upper[1] = e2;
  for (unsigned i=0; i<_ndims; ++i) 
    _length[i] = _upper[i]-_lower[i];
}

template<typename TYPE>
WxBox<TYPE>::WxBox(TYPE s1, TYPE e1, TYPE s2, TYPE e2, TYPE s3, TYPE e3)
  : _ndims(3)
{
  _lower[0] = s1; _upper[0] = e1;
  _lower[1] = s2; _upper[1] = e2;
  _lower[2] = s3; _upper[2] = e3;
  for (unsigned i=0; i<_ndims; ++i) 
    _length[i] = _upper[i]-_lower[i];
}

template<typename TYPE>
WxBox<TYPE>::WxBox(TYPE s1, TYPE e1, TYPE s2, TYPE e2, TYPE s3, TYPE e3, TYPE s4, TYPE e4)
  : _ndims(4)
{
  _lower[0] = s1; _upper[0] = e1;
  _lower[1] = s2; _upper[1] = e2;
  _lower[2] = s3; _upper[2] = e3;
  _lower[3] = s4; _upper[3] = e4;
  for (unsigned i=0; i<_ndims; ++i) 
    _length[i] = _upper[i]-_lower[i];
}

template<typename TYPE>
WxBox<TYPE>::WxBox(const WxBox& b)
{
  _ndims = b._ndims;
  for (unsigned i=0; i<_ndims; ++i)
  {
    _lower[i] = b._lower[i];
    _upper[i] = b._upper[i];
    _length[i] = b._length[i];
  }
}

template<typename TYPE>
WxBox<TYPE>&
WxBox<TYPE>::operator=(const WxBox& b) 
{
  if (this==&b)
    return *this;
    
  _ndims = b._ndims;
  for (unsigned i=0; i<_ndims; ++i)
  {
    _lower[i] = b._lower[i];
    _upper[i] = b._upper[i];
    _length[i] = b._length[i];
  }

  return *this;
}

template<typename TYPE>
WxBox<TYPE>::~WxBox()
{
}

template<typename TYPE>
void
WxBox<TYPE>::setup(const WxCryptSet& wxc)
{
  // extract stuff from crypt set

  std::vector<WxAny> lower = wxc.template get<std::vector<WxAny> >("Lower");
  std::vector<WxAny> upper = wxc.template get<std::vector<WxAny> >("Upper");

  for (unsigned i=0; i<_ndims; ++i)
  {
    _lower[i] = wx_any_cast<TYPE>(lower[i]);
    _upper[i] = wx_any_cast<TYPE>(upper[i]);
    _length[i] = _upper[i] - _lower[i];
  }
}

template<typename TYPE>
bool
WxBox<TYPE>::isEmpty() const
{
  for (unsigned i=0; i<_ndims; ++i)
    if (_length[i] <= 0)
      return true;
  return false;
}

template<typename TYPE>
TYPE
WxBox<TYPE>::area() const
{
  TYPE ar = (TYPE) 1;
  for (unsigned i=0; i<_ndims; ++i)
    ar *= _length[i];
  return ar;
}

template<typename TYPE>
bool
WxBox<TYPE>::operator==(const WxBox<TYPE>& b) const
{
  for (unsigned i=0; i<_ndims; ++i)
    if ((_lower[i] != b._lower[i]) || (_upper[i] != b._upper[i]))
      return false;
  return true;
}

template<typename TYPE>
WxBox<TYPE> 
WxBox<TYPE>::intersect(const WxBox<TYPE>& box) const {
  TYPE low[16], upp[16];
  for (unsigned i=0; i<_ndims; ++i) 
  {
    // lower is max of the two box's lower coordinates
    low[i] = this->lower(i) > box.lower(i) ? this->lower(i) : box.lower(i);
    // upper is min of the two box's upper coordinates
    upp[i] = this->upper(i) < box.upper(i) ? this->upper(i) : box.upper(i);
    // check if intersection is empty
    if (upp[i] <= low[i])
      // yes, so return empty box
      return WxBox<TYPE>(_ndims);
  }
  return WxBox<TYPE>(_ndims, low, upp);
}

template<typename TYPE>
WxBox<TYPE>
WxBox<TYPE>::intersectmultiblock(const WxBox<TYPE>& box, unsigned dir, int edge1, int edge2, bool recv, int npad) const {
  TYPE low[16], upp[16];

  int dirt, dirt1, dirt2;

  bool lowerBC = (edge1==0)&&(edge2==1);
  bool upperBC = (edge1==1)&&(edge2==0);
  bool boxnotamatch, thisnotamatch, notamatch;

  for (unsigned i=0; i<_ndims; ++i)
  {
	  if (i!=dir)
	  {
		  // lower is max of the two box's lower coordinates
		  low[i] = this->lower(i) > box.lower(i) ? this->lower(i) : box.lower(i);
		  // upper is min of the two box's upper coordinates
		  upp[i] = this->upper(i) < box.upper(i) ? this->upper(i) : box.upper(i);
	  }
	  else
	  {
		  boxnotamatch = (box.lower(i)==box.upper(i));
		  thisnotamatch = (this->lower(i)==this->upper(i));
		  notamatch = boxnotamatch || thisnotamatch;
//		  std::cout << "the box is not a match = "<< boxnotamatch << std::endl;
//		  std::cout << "the this is not a match = "<< thisnotamatch << std::endl;
//		  std::cout << "not a match = "<< notamatch << std::endl;

		  if (notamatch)
		  {
			  low[i] = 0;
			  upp[i] = 0;
//			  std::cout << "not a match executing" << std::endl;
		  }
		  else
		  {
			  if (not(recv))
			  {
				  low[i] = this->lower(i);
				  upp[i] = this->upper(i);
			  }
			  else
			  {
				  if (upperBC)
				  {
					  low[i] = this->lower(i)+npad;
					  upp[i] = this->upper(i)+npad;
				  }
				  if (lowerBC)
				  {
					  low[i] = this->lower(i)-npad;
					  upp[i] = this->upper(i)-npad;
				  }
			  }
		  }
	  }
	  // check if intersection is empty
//	  std::cout << "(_ndims, low["<<i<<"], upp["<<i<<"]) = (" <<_ndims << "," << low[i] <<"," << upp[i] << ")" << std::endl;

	  if (upp[i] <= low[i])
	  {
		  // yes, so return empty box
//		  std::cout << "intersectmultiblock returns an empty box because upp["<< i << "] = "<< upp[i] << "is equal to or less than low["<< i << "] = "<< low[i] <<"." << std::endl;
		  return WxBox<TYPE>(_ndims);
	  }
//	  std::cout<< "intersectmultiblock returns the following:  (_ndims, low["<<i<<"], upp["<<i<<"]) = (" <<_ndims << "," << low[i] <<"," << upp[i] << ")" << std::endl;
  }
  return WxBox<TYPE>(_ndims, low, upp);
}

template<typename TYPE>
WxBox<TYPE>
WxBox<TYPE>::extend(const TYPE low[], const TYPE upp[]) const
{
  TYPE newLow[16], newUpp[16];
  for (unsigned i=0; i<_ndims; ++i)
  {
    newLow[i] = _lower[i] - low[i];
    newUpp[i] = _upper[i] + upp[i];
  }
  return WxBox<TYPE>(_ndims, newLow, newUpp);
}

template<typename TYPE>
WxBox<TYPE>
WxBox<TYPE>::extDim(TYPE low, TYPE upp) const
{
  TYPE newLow[16], newUpp[16];
  newLow[0] = low; newUpp[0] = upp;
  for (unsigned i=0; i<_ndims; ++i)
  {
    newLow[i+1] = _lower[i];
    newUpp[i+1] = _upper[i];
  }
  return WxBox<TYPE>(_ndims+1, newLow, newUpp);
}

template<typename TYPE>
bool
operator==(const WxBox<TYPE>& ra, const WxBox<TYPE>& rb)
{
  if (ra.ndims() != rb.ndims())
    return false;
  for (unsigned i=0; i<ra.ndims(); ++i)
    if ( (ra.lower(i) != rb.lower(i)) || (ra.upper(i) != rb.upper(i)) )
      return false;
  return true;
}

template<typename TYPE>
bool
WxBox<TYPE>::contains(TYPE coord[]) const
{
  for (unsigned i = 0; i < _ndims; ++i)
    if ( (coord[i] < _lower[i]) || (coord[i] >= _upper[i]) )
      return false;
  return true;
}

// instantiations
template class WxBox<int>;
//template class WxBox<float>;
template class WxBox<double>;
