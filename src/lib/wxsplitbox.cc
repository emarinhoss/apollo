// WarpX lib includes
#include "wxsplitbox.h"
#include "wxrange.h"
#include "wxsequencer.h"
#include "wxexcept.h"
#include "wxlogger.h"
#include "wxlogstream.h"

// std includes
#include <cmath>
#include <iostream>

template <typename TYPE>
WxSplitBox<TYPE>::WxSplitBox(WxMsgBase *msgBase, const WxBox<TYPE>& box, bool pdirs[])
  : NDIMS(box.ndims()), 
    _useCount(new int),
    _data(new WxSplitBoxData(msgBase->numProcs(), box, msgBase))
{
  *_useCount = 1;
  if (pdirs)
  {
    for (unsigned i=0; i<box.ndims(); ++i)
      _data->periods[i] = pdirs[i];
  }
  else
  {
    for (unsigned i=0; i<box.ndims(); ++i)
      _data->periods[i] = false;
  }
  _divideBox(msgBase->numProcs(), box);
  _data->nboxes = msgBase->numProcs();
  _duplicateBoxes();
}

template <typename TYPE>
WxSplitBox<TYPE>::WxSplitBox(WxMsgBase *msgBase, unsigned subBox[], const WxBox<TYPE>& box, bool pdirs[])
  : NDIMS(box.ndims()), _useCount(new int), _data(new WxSplitBoxData(1, box, msgBase))
{
  *_useCount = 1;
  if (pdirs)
  {
    for (unsigned i=0; i<box.ndims(); ++i)
      _data->periods[i] = pdirs[i];
  }
  else
  {
    for (unsigned i=0; i<box.ndims(); ++i)
      _data->periods[i] = false;
  }

  std::vector<WxBox<TYPE> > boxes = _breakBoxes(subBox, box);
  // insert boxes into map
  unsigned count = 0;
  typename std::vector<WxBox<TYPE> >::const_iterator i;
  for (i = boxes.begin(); i != boxes.end(); ++i)
    _data->boxMap.insert( BoxPair_t(count++, *i) );
  _data->nboxes = count;
  _duplicateBoxes();

}

template <typename TYPE>
WxSplitBox<TYPE>::~WxSplitBox() 
{
  if (--*_useCount == 0)
  {
    delete _data;
    delete _useCount;
  }
}

template <typename TYPE>
WxSplitBox<TYPE>::WxSplitBox(const WxSplitBox<TYPE>& bs)
  : NDIMS(bs.NDIMS), _data(bs._data)
{
  ++*bs._useCount;
  _useCount = bs._useCount;
}

template <typename TYPE>
WxSplitBox<TYPE>&
WxSplitBox<TYPE>::operator=(const WxSplitBox<TYPE>& bs)
{
  NDIMS = bs.NDIMS;
  ++*bs._useCount;
  if (--*_useCount == 0)
  {
    delete _data;
  }
  _data = bs._data;
  _useCount = bs._useCount;
  return *this;
}

template<typename TYPE>
WxBox<TYPE> 
WxSplitBox<TYPE>::getBox(unsigned n) const 
{
  typename BoxMap_t::const_iterator i = _data->boxMap.find(n);
  if (i != _data->boxMap.end())
    return (*i).second;
  // not found: throw an exception
  throw WxExcept("Box not found");
}

template<typename TYPE>
unsigned
WxSplitBox<TYPE>::getBoxRank(unsigned n) const {
  typename BoxRank_t::const_iterator i = _data->boxRank.find(n);
  if (i != _data->boxRank.end())
    return (*i).second;
  // not found: throw an exception
  throw WxExcept("Box not found");
}

template <typename TYPE>
void
WxSplitBox<TYPE>::_divideBox(unsigned num, const WxBox<TYPE>& decompBox) 
{
  //  begin by determining the number of breaks in the first N-1 dimension
  double side = pow((double) num, (double) 1 / NDIMS);
  unsigned sbounds[2];
  sbounds[0] = (unsigned) std::floor(side);
  sbounds[1] = (unsigned) std::ceil(side);
  unsigned boundtest = 1;
  unsigned breaks[16];
  for (unsigned i = 0; i < NDIMS; ++i) 
  {
    boundtest *= sbounds[0];
  }
  //  This does a series of comparisons to figure out how to distribute boxes
  //  in the first N-1 dimensions
  for (unsigned i = 0; i < NDIMS; ++i) 
  {
    boundtest = boundtest * sbounds[1] / sbounds[0];
    if (num <= boundtest) {
      for (unsigned j = 0; j < i; ++j) 
      {
        breaks[j] = sbounds[1];
      }
      for (unsigned j = i; j < NDIMS - 1; ++j) 
      {
        breaks[j] = sbounds[0];
      }
      break;
    }
  }
  breaks[NDIMS-1] = 1;
  //  Now break down the box into a series of strips
  std::vector< WxBox<TYPE> > intermediate = _breakBoxes(breaks, decompBox);
  //  Do the final breakdown, the first rmdr boxes have sbrk + 1
  //  breaks in it, while the rest have sbrk breaks in it.  rmdr is
  //  the remainder of num strips into num processes.  sbrk is the
  //  quotient of num strips into num processes.  Also rebalance boxes.
  unsigned rmdr = num % intermediate.size();
  unsigned sbrk = num / intermediate.size();
  unsigned dims1[16], *dimsary;
  for (unsigned i = 0; i < NDIMS - 1; ++i) 
  {
    dims1[i] = breaks[i];
  }
  //  This whole godawful messy thing rebalances boxes to all have
  //  ~ equal areas.
  dimsary = new unsigned[intermediate.size()];
  for (unsigned i = 0; i < intermediate.size(); ++i) 
  {
    dimsary[i] = i < rmdr ? sbrk + 1 : sbrk;
  }
  for (unsigned i = 0; i < NDIMS - 1; ++i) 
  {
    unsigned step = 1;
    for (unsigned j = 0; j < NDIMS - 1 - i; ++j) 
    {
      step *= dims1[j];
    }
    unsigned substep = 1;
    for (unsigned j = 0; j < NDIMS - 2 - i; ++j) 
    {
      substep *= dims1[j];
    }
    for (unsigned n = 0; n < intermediate.size(); n += step) 
    {
      unsigned weightsum = 0;
      for (unsigned j = n; j < n + step; ++j) 
      {
        weightsum += dimsary[j];
      }
      unsigned rollingsum = 0;
      for (unsigned j = n; j < n + step; j += substep) 
      {
        unsigned stepsum = 0;
        for (unsigned k = j; k < j + substep; ++k) 
        {
          stepsum += dimsary[k];
        }
        rollingsum += stepsum;
        TYPE newUpper = rollingsum * decompBox.length(NDIMS-2-i) / weightsum;
        for (unsigned k = j; k < j + substep; ++k) 
        {
          intermediate[k].upper(NDIMS-2-i, newUpper);
          if (j + substep < n + step) 
          {
            intermediate[k+substep].lower(NDIMS-2-i, newUpper);
          }
        }
      }
    }
  }
  delete[] dimsary;
  for (unsigned i = 0; i < NDIMS - 1; ++i) 
  {
    breaks[i] = 1;
  }
  unsigned count = 0;
  for (unsigned i = 0; i < intermediate.size(); ++i) 
  {
    breaks[NDIMS-1] = i < rmdr ? sbrk + 1 : sbrk;
    std::vector< WxBox<TYPE> > temp;
    temp = _breakBoxes(breaks, intermediate[i]);
    // insert boxes
    typename std::vector<WxBox<TYPE> >::const_iterator j;
    for (j = temp.begin(); j != temp.end(); ++j)
      _data->boxMap.insert( BoxPair_t(count++, *j) );
  }
}

template <typename TYPE>
void
WxSplitBox<TYPE>::_duplicateBoxes() 
{
  BoxMap_t realMap(_data->boxMap.begin(), _data->boxMap.end());
  typename BoxMap_t::const_iterator bi;
  for (bi = realMap.begin(); bi != realMap.end(); ++bi)
  {
    BoxRPair_t insertme((*bi).first, (*bi).first);
    _data->boxRank.insert(insertme);
  }
  int itrMin[16], itrMax[16];
  for (unsigned i = 0; i < NDIMS; ++i) 
  {
    itrMin[i] = _data->periods[i] ? -1 : 0;
    itrMax[i] = _data->periods[i] ? 2 : 1;
  }
  WxBox<int> shift(NDIMS, itrMin, itrMax);
  WxSequencer<> seq(shift);
  unsigned nTotBoxes = _data->nboxes;
  while ( seq.step() ) {
    //  Skip the case of zero offset
    bool skipStep = true;
    for (unsigned i = 0; i < NDIMS; ++i) 
    {
      if (seq.indices()[i] != 0) 
      {
        skipStep = false;
        break;
      }
    }
    //  Shift boxes, renumber, and add to the box map
    if (!skipStep) 
    {
      TYPE lext[16], uext[16], dist[16];
      for (unsigned i = 0; i < NDIMS; ++i) 
      {
        dist[i] = _data->box.length(i);
        lext[i] = -1 * seq.indices()[i] * dist[i];
        uext[i] =  1 * seq.indices()[i] * dist[i];
      }
      typename BoxMap_t::const_iterator bi2;
      for (bi2 = realMap.begin(); bi2 != realMap.end(); ++bi2) 
      {
        BoxPair_t insertme(nTotBoxes, (*bi2).second.extend(lext, uext));
        BoxRPair_t insertmetoo(nTotBoxes, (*bi2).first);
        _data->boxMap.insert(insertme);
        _data->boxRank.insert(insertmetoo);
        nTotBoxes++;
      }
    }
  }
}

template <typename TYPE>
std::vector<WxBox<TYPE> >
WxSplitBox<TYPE>::_breakBoxes(unsigned subBox[16], const WxBox<TYPE>& box) 
{
  std::vector<WxBox<TYPE> > boxes;
  // compute total number of boxes
  unsigned nbrks = 1;
  for (unsigned i=0; i<NDIMS; ++i)
    nbrks *= subBox[i];
  // construct boxes
  unsigned itr[16];
  for (unsigned i=0; i<NDIMS; ++i)
    itr[i] = 0;
  TYPE newLower[16], newUpper[16];
  for (unsigned i=0; i<nbrks; ++i)
  {
    for (unsigned j=0; j<NDIMS; ++j) 
    {
      newLower[j] = box.lower(j) + itr[j] * box.length(j) / subBox[j];
      newUpper[j] = box.lower(j) + (itr[j] + 1) * box.length(j) / subBox[j];
    }
    boxes.push_back(WxBox<TYPE>(NDIMS, newLower, newUpper));

    for (unsigned j=0; j<NDIMS; ++j) 
    {
      if (++itr[j] >= subBox[j])
        itr[j] = 0;
      else
        break;
    }
  }
  return boxes;
}

template <typename TYPE>
unsigned
WxSplitBox<TYPE>::getRank(TYPE coords[]) const 
{
  typename BoxMap_t::const_iterator i, iend;
  iend = _data->boxMap.end();
  for (i = _data->boxMap.begin(); i != iend; ++i) {
    if ((*i).second.contains(coords)) {
      return (*i).first;
    }
  }
  //  Doesn't exist in the decomp: throw a tde.
  WxExcept wxe("WxSplitBox::getRank: ");
  wxe << "Point " << coords[0];
  for (unsigned i = 1; i < NDIMS; ++i)
    wxe << "," << coords[i];
  wxe << " is not in decomposition.";
  throw wxe;
}

// instantiations
template class WxSplitBox<int>;
template class WxSplitBox<float>;
template class WxSplitBox<double>;
