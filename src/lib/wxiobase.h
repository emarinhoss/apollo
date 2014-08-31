#ifndef __wxiobase__
#define __wxiobase__

// WarpX includes
#include "wxiotmpl.h"
#include "wxtypelist.h"
#include "wxdatatypes.h"
#include "wxexcept.h"

// std includes
#include <vector>

/**
 * Provides an abstract interface for access to hierachical datasets
 */
class WxIoBase {

  public:

/**
 * Virtual destructor
 */
    virtual ~WxIoBase();

/**
 * Set the base name
 *
 * @param bn the base name
 */
    void setBaseName(const std::string& bn) {
      baseName = bn;
    }

/**
 * Get the base name
 *
 * @return the base name.  Used to compose the file name.
 */
    std::string getBaseName() {
      return baseName;
    }

/**
 * Set the dump number
 *
 * @param d the dump number.  Not used if negative.
 */
    void setDumpNo(int d) {
      dumpNo = d;
    }

/**
 * Get the dump number
 *
 * @return the dump number
 */
    int getDumpNo() {
      return dumpNo;
    }

/**
 * Create a file.
 *
 * @param fileName the name for the file.  Assumed to be rw.
 * @return node for the file
 */
    virtual WxIoNodeType createFile(const std::string& fileName) = 0;

/**
 * Create a file with base name prepended and dump and suffix appended.
 *
 * @param dataName the name for the file.  Assumed to be rw.
 * @return node for the file
 */
    virtual WxIoNodeType createDumpFile(const std::string& dataName);

/**
 * Open a file.
 *
 * @param fileName the name for the file
 * @param perms: the read and write permissions.  "r" or "rw"
 * @return node for the file
 */
    virtual WxIoNodeType openFile(const std::string& fileName, const std::string& perms) = 0;

/**
 * Open a file with the base name prepended and the dump and suffix
 * appended.
 *
 * @param dataName used to compose the name of the file
 * @param perms: the read and write permissions.  "r" or "rw"
 * @return node for the file
 */
    virtual WxIoNodeType openDumpFile(const std::string& dataName,
	const std::string& perms);

/**
 * Close a file node
 */
    virtual void closeFile(WxIoNodeType fileNode) = 0;

/**
 * Get the dump name for a given data name
 */
    virtual std::string getDumpFileName(const std::string& dataName);

/**
 * Create an empty group
 *
 * @param node the node to write under
 * @param dataName the name of the data
 *
 * @return the node to the written data.  This is not closed.
 */
    virtual WxIoNodeType createGroup(WxIoNodeType node, const std::string& dataName) const = 0;

/**
 * Open a group
 *
 * @param node the node to look in
 * @param dataName the name of the data to open
 *
 * @return the node to the read data
 */
    virtual WxIoNodeType openGroup(WxIoNodeType node, const std::string& dataName) const = 0;

/**
 * Create an empty node
 *
 * @param node the node to write under
 * @param dataName the name of the data
 *
 * @return the node to the written data.  This is not closed.
 */
    virtual WxIoNodeType createDataSet(WxIoNodeType node, const std::string& dataName) const = 0;

/**
 * Open a node
 *
 * @param node the node to look in
 * @param dataName the name of the node to open
 *
 * @return the node to the read data
 */
    virtual WxIoNodeType openDataSet(WxIoNodeType node, const std::string& dataName) const = 0;

/**
 * Write a new data set under a node.
 *
 * @param node the node to write under
 * @param dataName the name of the data
 * @param dataSetSize vector of the sizes of the entire data set
 * @param dataSetBeg first index for each direction
 * @param dataSetLen length of data to be written for each direction
 * @param data the data to be written
 * @return the node to the written data.  The node is not closed.
 */
    template <class DATATYPE> WxIoNodeType writeDataSet(WxIoNodeType node,
	const std::string& dataName, const std::vector<size_t>& dataSetSize,
	const std::vector<size_t>& dataSetBeg,
	const std::vector<size_t>& dataSetLen, const DATATYPE* data) {
      return getIoPtr<DATATYPE>()->writeDataSet(node, dataName, dataSetSize,
	dataSetBeg, dataSetLen, data);
    }

/**
 * Read a data set under a node.
 *
 * @param node the node to write under
 * @param dataName the name of the data
 * @param dataSetBeg first index for each direction
 * @param dataSetLen length of data to be written for each direction
 * @param data the data to be read.
 * @return the node to the read data. The node is not closed.
 */
    template <class DATATYPE> WxIoNodeType readDataSet(WxIoNodeType node,
	const std::string& dataName, const std::vector<size_t>& dataSetBeg,
	const std::vector<size_t>& dataSetLen, DATATYPE* data) {
      return getIoPtr<DATATYPE>()->readDataSet(node, dataName, dataSetBeg,
	dataSetLen, data);
    }

/**
 * Close a data set.
 *
 * @param node the node to be closed
 */
    virtual void closeDataSet(WxIoNodeType node) const = 0;

/**
 * Write an attribute.
 *
 * @param node the node to write under
 * @param attribName the name of the attribute
 * @param attrib the attribute to be written
 */
    template <class DATATYPE> void writeAttribute(WxIoNodeType node,
	const std::string& attribName, const DATATYPE& attrib) {
        getIoPtr<DATATYPE>()->writeAttribute(node, attribName, attrib);
    }

/**
 * Write a vector of attributes.
 *
 * @param node the node to write under
 * @param attribName the name of the attribute
 * @param attrib the vector attribute to be written
 */
    template <class DATATYPE> void writeVecAttribute(WxIoNodeType node,
	const std::string& attribName, const std::vector<DATATYPE>& attrib) {
        getIoPtr<DATATYPE>()->writeVecAttribute(node, attribName, attrib);
    }

/**
 * Write a string attribute
 *
 * @param node the node to write under
 * @param attribName the name of the attribute
 * @param attrib the string value of the attribute
 */
    virtual void writeStrAttribute(WxIoNodeType node, const std::string& attribName, const std::string& attrib) const = 0;

/**
 * Read an attribute.
 *
 * @param node the node to write under
 * @param attribName the name of the attribute
 * @param attrib the attribute to be read
 */
    template <class DATATYPE> void readAttribute(WxIoNodeType node,
	const std::string& attribName, DATATYPE& attrib) {
        getIoPtr<DATATYPE>()->readAttribute(node, attribName, attrib);
    }

/**
 * Read an attribute.
 *
 * @param node the node to write under
 * @param attribName the name of the attribute
 * @param attrib the vector attribute to be read
 */
    template <class DATATYPE> void readVecAttribute(WxIoNodeType node,
	const std::string& attribName, std::vector<DATATYPE>& attrib) {
        getIoPtr<DATATYPE>()->readVecAttribute(node, attribName, attrib);
    }

  protected:

/**
 * Default constructor - no dump, basename or suffix.
 *
 */
    WxIoBase();

/**
 * Constructor is protected, as this class cannot be made standalone.
 * Construct with suffix only.
 *
 * @param sfx the suffix for file names
 */
    WxIoBase(const std::string& sfx);

/**
 * Constructor is protected, as this class cannot be made standalone.
 *
 * @param bn the base for file names
 * @param d the dump number for file names
 * @param sfx the suffix for file names
 */
    WxIoBase(const std::string& bn, int d, const std::string& sfx);

/**
 * Add a file to the list of open files.
 *
 * @param node the file node
 */
    virtual void addOpenFile(WxIoNodeType node);

/**
 * Remove a file from the list of open files.
 *
 * @param node the file node
 */
    virtual void removeOpenFile(WxIoNodeType node);

/**
 * Close any files currently open.
 */
    virtual void closeOpenFiles();

/**
 * Add a new io object. The derived class should call this to setup
 * WxIoBase properly.
 */
    template <typename T>
    void addIo(const WxIoTmpl<T>* b) {
      wxTypeMapExtract<T>(ioTypeMap).ioPtr = b;
    }

  private:

/** Private copy constructor to prevent use */
    WxIoBase(const WxIoBase&);

/** Private assignment to prevent use */
    WxIoBase& operator=(const WxIoBase&);

/**
 * Get an io object of the needed type
 */
    template <typename T>
    const WxIoTmpl<T>* getIoPtr() {
        const WxIoTmpl<T>* r = wxTypeMapExtract<T>(ioTypeMap).ioPtr;
      if (r) return r;
      WxExcept wxe;
      wxe << "I/O type not set properly";
      throw wxe;
    }

    std::vector<WxIoNodeType> openFiles; // List of all open files
    std::string baseName; // base name for output
    int dumpNo; // dump number
    std::string suffix; // suffix for output

//  Note: this struct is public as xlC doesn't like internal
//  private structs
  public:

    // container class for all io classes
    template <typename T>
    struct WxIoContainer {
        WxIoContainer() : ioPtr(0) {
        }
        virtual ~WxIoContainer() {
            delete ioPtr;
        }
        // Points to a derived class of WxIoTmpl<T> 
        const WxIoTmpl<T>* ioPtr;
    };

    // Objects of type WxIoTypeMap_t inherit from all WxIoTmpl<T>
    // where T belongs to the WxIoTypelist_t. Thus it acts like a
    // container for all io objects in the system.
    typedef WxTypeMap<WxDataTypes_t, WxIoContainer> WxIoTypeMap;

    WxIoTypeMap ioTypeMap; // Container of io
};

#endif // __wxiobase__

