#ifndef __wxhdf5io__
#define __wxhdf5io__

// WarpX lib includes
#include "wxiobase.h"

#ifdef _DO_USE_MPI_
#include <mpi.h>
#else
/** Communicator set to void if not using MPI */
#define MPI_Comm void*
/** Info set to void if not using MPI */
#define MPI_Info void*
#endif

/**
 * WxHdf5Io is the interface for the HDF5 implementation of HDF5.
 */
class WxHdf5Io : public WxIoBase {

  public:

/**
 * Constructor creates the individual templated writers
 *
 * @param mc the MPI communicator
 * @param mi info for the MPI communicator
 */
    WxHdf5Io(MPI_Comm mc, MPI_Info mi);

/**
 * Constructor creates the individual templated writers
 *
 * @param bn base name for the dump
 * @param d the dump number
 * @param mc the MPI communicator
 * @param mi info for the MPI communicator
 */
    WxHdf5Io(const std::string& bn, int d, MPI_Comm mc, MPI_Info mi);

/**
 * Virtual destructor
 */
    virtual ~WxHdf5Io();

/**
 * Create a file.
 *
 * @param fileName the name for the file.  Assumed to be rw.
 *
 * @return node for the file
 */
    virtual WxIoNodeType createFile(const std::string& fileName);

/**
 * Open a file.
 *
 * @param fileName the name for the file
 * @param perms: the read and write permissions.  "r" or "rw"
 *
 * @return node for the file
 */
    virtual WxIoNodeType openFile(const std::string& fileName,
        const std::string& perms);

/**
 * Create an empty group
 *
 * @param node the node to write under
 * @param dataName the name of the data
 *
 * @return the node to the written data
 */
    virtual WxIoNodeType createGroup(WxIoNodeType node, const std::string& dataName) const;

/**
 * Open a group
 *
 * @param node the node to look in
 * @param dataName the name of the data to open
 *
 * @return the node to the read data
 */
    virtual WxIoNodeType openGroup(WxIoNodeType node, const std::string& dataName) const;

/**
 * Create an empty node
 *
 * @param node the node to write under
 * @param dataName the name of the data
 *
 * @return the node to the written data
 */
    virtual WxIoNodeType createDataSet(WxIoNodeType node, const std::string& dataName) const;

/**
 * Open a node
 *
 * @param node the node to look in
 * @param dataName the name of the node to open
 *
 * @return the node to the read data
 */
    virtual WxIoNodeType openDataSet(WxIoNodeType node, const std::string& dataName) const;

/**
 * Get the top node = file node
 */
    virtual void closeFile(WxIoNodeType fileNode);

/**
 * Close a data set.
 *
 * @param node the node to be closed
 */
    virtual void closeDataSet(WxIoNodeType node) const;

/**
 * Write a string attribute
 *
 * @param node the node to write under
 * @param attribName the name of the attribute
 * @param attrib the string value of the attribute
 */
    void writeStrAttribute(WxIoNodeType node, const std::string& attribName, const std::string& attrib) const;

  protected:

  private:

/** Private copy constructor to prevent use */
    WxHdf5Io(const WxHdf5Io&);

/** Private assignment to prevent use */
    WxHdf5Io& operator=(const WxHdf5Io&);

/**
 * Set up the templated IO objects and the mpi stuff
 */
    void setup(MPI_Comm mc, MPI_Info mi);

/** The communicator */
    MPI_Comm mpiComm;

/** The mpi info */
    MPI_Info mpiInfo;

};

#endif // __wxhdf5io__

