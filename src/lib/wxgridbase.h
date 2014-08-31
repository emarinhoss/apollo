#ifndef __wxgridbase__
#define __wxgridbase__

#include <wxobject.h>

class WxGridBase : public WxObject
{
  public:
/**
 * Setup grid using supplied crypset.
 *
 * @param wxc Cryptset using which the object is set up.
 */
    virtual void setup(const WxCryptSet& wxc);

/**
 * Initialize the grid.
 */
    virtual void init();

/**
 * Load grid from file.
 *
 * @param grpNode group node to read from
 */
    virtual void load(const WxIoNodeType& grpNode);

/**
 * Dump grid to file.
 *
 * @param io I/O object to use for writing
 * @param grpNode group node to write to
 */
    virtual void dump(WxIoBase& io, WxIoNodeType& grpNode);

  private:
    
};

#endif // __wxgridbase__
