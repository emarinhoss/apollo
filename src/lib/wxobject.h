#ifndef __wxobject__h__
#define __wxobject__h__

// WarpX lib includes
#include <wxiobase.h>
#include <wxmsgbase.h>
#include <wxcryptset.h>

// PETSc includes
#include <petsc.h>

// std includes
#include <iostream>
#include <string>

/**
 * WxObject is a base class for WarpX classes which need to go through
 * a creation/destruction cycle within the simulation. Objects derived
 * from this class can provide several methods which are used to
 * create and destroy them. 
 *
 * The objects in WarpX are created in a three-step process. First,
 * the setup() method is called which should read in data from the
 * cryptset. Second, either the init() or the load() method is
 * called. The init() method is not called if the object is being
 * restored from an output file. Instead, the load() method is
 * called. Finally, the finishBuild() method is called.
 * 
 * The init() method should initialize the object using the data read
 * in by the setup() method. Often, for simple objects, this method is
 * not needed and all the work can be done in the setup() method.
 *
 * The load() method should read data from the HDF5 file and restore
 * itself. This method is passed a pointer to the HDF5 group from
 * which the object should read data to restore itself.
 *
 * Finally, the finishBuild() method is called which should complete
 * the initialization of the object. After this call is complete the
 * object should be ready for use.
 *
 * Simple objects may not need init(), load(), dump() and
 * finishBuild() methods. In this case they need not be provided by
 * the derived classes.
 *
 * Objects which provide a load() need to provide a dump()
 * method. This method should write its data out to the HDF5 group
 * provided. It should write enough data so that it can reconstruct
 * itself when the load() method is called.
 */
class WxObject
{
  public:
/** Create object */
    WxObject();

/** 
 * Create object with given name
 * 
 * @param name Name of object
 */
    WxObject(const std::string& name);
    
/** Dtor: destroy object */
    virtual ~WxObject();

/**
 * Set the I/O pointer for use in object
 *
 * @param io I/O object to use
 */
    void setIo(WxIoBase& io);

/**
 * Set the msg pointer for use in object
 *
 * @param msg object to use
 */
    void setMsg(WxMsgBase& msg);

/**
 * Return reference to I/O object
 *
 * @return reference to I/O object
 */
    WxIoBase& getIo();

/**
 * Return reference to msg object
 *
 * @return reference to msg object
 */
    WxMsgBase& getMsg();

/**
 * Setup object using supplied crypset. This method should read in
 * parameters needed to initialize it.
 *
 * @param wxc Cryptset using which the object is set up.
 */
    virtual void setup(const WxCryptSet& wxc);

/**
 * Initialize the object. This is called after the setup() method is
 * called. This method is not called if the object is being restored
 * from an output file. In this case the load() method is called
 * instead.
 */
    virtual void init(PetscReal dt, Vec X);

/**
 * Finish building the object. This is called after the init() or
 * load() methods have been called. In this the object should finish
 * building itself and be ready for use.
 */
    virtual void finishBuild();

/**
 * Load object from file. This method is passed a group node from
 * which its data should be read and then constructed.
 *
 * @param io I/O object to use for writing
 * @param grpNode group node to read from
 */
    virtual void load(WxIoBase& io, const WxIoNodeType& grpNode);

/**
 * Dump object to file. This methods is passed a group node to which
 * is should write its data to. Enough data should be written so that
 * it can reconstruct itself using the load() method.
 *
 * @param io I/O object to use for writing
 * @param grpNode group node to write to
 */
    virtual void dump(WxIoBase& io, WxIoNodeType& grpNode);
 
/**
 * Get name of object.
 *
 * @return Name of object
 */
    std::string name() const;

/**
 * Set object's name
 */
    void setName(const std::string& nm);

  private:
/** Name of object */
    std::string _name;
/** Pointer to I/O object */
    WxIoBase *_io;
/** Pointer to messaging object */
    WxMsgBase *_msg;
};

#endif //  __wxobject__h__
