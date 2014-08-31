#ifndef __wxiotmpl__
#define __wxiotmpl__

// std includes
#include <string>
#include <vector>

/**
 * Provides a means for derived messengers to return implimentation
 * specific message status flags and data. The returned type is an
 * opaque pointer and should not be directly fiddled around with.
 */
struct WxIoNodeTypev {

/** Destructor */
    virtual ~WxIoNodeTypev() {}

/**
 * A comparison operator
 *
 * @param v the node to which we are comparing ourself
 * @return a boolean specifying whether or not we are equal
 */
    virtual bool operator==(const WxIoNodeTypev& v) = 0;

};

/**
 * A simple typedef to easily refer to an WxIoNodeType pointer w/o
 * having to know about it
 */
typedef WxIoNodeTypev* WxIoNodeType;

/**
 * WxIoTmpl is the base class for access to a hierarchical file
 * system with groups, data sets, and attributes for those datasets.
 * The exemplar is HDF5, but one may eventually other systems, like
 * netCDF, PDB, or ?
 *
 * This is the base class.  It will provide for output of the
 * basic data structures of the STL library and boost objects.
 * Derived classes will allow for output of Facets objects, such
 * as Facets arrays, whether distributed or not.
 */
template <class DATATYPE>
class WxIoTmpl {

  public:

/**
 * Virtual destructor
 */
    virtual ~WxIoTmpl() {
    }

/**
 * Write a new data set under a node
 *
 * @param node the node to write under
 * @param dataName the name of the data
 * @param dataSetSize vector of the sizes of the entire data set
 * @param dataSetBeg first index for each direction
 * @param dataSetLen length of data to be written for each direction
 * @param data the data to be written
 *
 * @return the node to the written data
 */
    virtual WxIoNodeType writeDataSet(WxIoNodeType node,
        const std::string& dataName, const std::vector<size_t>& dataSetSize,
        const std::vector<size_t>& dataSetBeg,
        const std::vector<size_t>& dataSetLen, const DATATYPE* data) const = 0;

/**
 * Read a new data set under a node
 *
 * @param node the node to write under
 * @param dataName the name of the data
 * @param dataSetBeg first index for each direction
 * @param dataSetLen length of data to be written for each direction
 * @param data the data to be written
 *
 * @return the node to the written data
 */
    virtual WxIoNodeType readDataSet(WxIoNodeType node,
        const std::string& dataName, const std::vector<size_t>& dataSetBeg,
        const std::vector<size_t>& dataSetLen, DATATYPE* data) const = 0;

/**
 * Write an attribute.
 *
 * @param node the node to write under
 * @param attribName the name of the attribute
 * @param attrib the attribute to be written
 */
    virtual void writeAttribute(WxIoNodeType node,
        const std::string& attribName, const DATATYPE& attrib) const = 0;

/**
 * Write an attribute.
 *
 * @param node the node to write under
 * @param attribName the name of the attribute
 * @param attrib the vector attribute to be written
 */
    virtual void writeVecAttribute(WxIoNodeType node,
        const std::string& attribName, const std::vector<DATATYPE>& attrib) const = 0;

/**
 * Read an attribute.
 *
 * @param node the node to which this attribute belongs
 * @param attribName the name of the attribute
 * @param attrib the attribute to be read
 */
    virtual void readAttribute(WxIoNodeType node,
        const std::string& attribName, DATATYPE& attrib) const = 0;

/**
 * Read an attribute.
 *
 * @param node the node to which this attribute belongs
 * @param attribName the name of the attribute
 * @param attrib the vector attribute to be read
 */
    virtual void readVecAttribute(WxIoNodeType node,
        const std::string& attribName, std::vector<DATATYPE>& attrib) const = 0;

  protected:

/**
 * Constructor is protected, as this class cannot be made standalone.
 */
    WxIoTmpl() {
    }

  private:

/** Private copy constructor to prevent use */
    WxIoTmpl(const WxIoTmpl<DATATYPE>&);

/** Private assignment to prevent use */
    WxIoTmpl<DATATYPE>& operator=(const WxIoTmpl<DATATYPE>&);

};

#endif // __wxiotmpl__

