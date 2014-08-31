#ifndef __wxhdf5iotmpl__
#define __wxhdf5iotmpl__

// WarpX lib includes
#include "wxiotmpl.h"

// HDF5 includes
#define H5_USE_16_API
#include <hdf5.h>

// std includes
#include <vector>

/**
 * HDF5 specific i/o node wrapper
 */
struct WxHdf5NodeTypev : public WxIoNodeTypev {

/**
 * The constructor
 *
 * @param node the node number of this node
 * @param parent the parent of this node
 */
    WxHdf5NodeTypev(long node=0, WxHdf5NodeTypev* parent=0)
            : node(node), parent(parent), nodeType(1) {
        if (parent)
            parent->children.push_back(this);
    }

/** Destructor */
    virtual ~WxHdf5NodeTypev() {
        // loop over children killing them
        unsigned n = children.size();
        for (unsigned i=0; i<n; ++i) {
            WxIoNodeTypev *ptr = children.back();
            delete ptr;
        }
        if (parent) {
            if (nodeType == 1) {
                H5Dclose(node);
            } else if (nodeType == 0) {
                H5Gclose(node);
            }
            // remove us from parent
            std::vector<WxHdf5NodeTypev*>::iterator itr;
            for (itr = parent->children.begin(); itr != parent->children.end(); ++itr) {
                if (**itr == *this) {
                    parent->children.erase(itr);
                    break;
                }
            }
        }
        else {
            // we are file node
            H5Fclose(node);
        }
    }

/**
 * Sets the type of this node
 *
 * @param the type of this node
 */
    void setNodeType(unsigned type) {
        nodeType = type;
    }

/**
 * Comparison operator
 *
 * @param v the node to which we are comparing ourself
 * @return boolean specifying whether we are the same
 */
    bool operator==(const WxIoNodeTypev& v) {
        return node == static_cast<const WxHdf5NodeTypev&>(v).node;
    }

/** This node's node number */
    long node;
/** This node's parent node */
    WxHdf5NodeTypev *parent;
/** This node's list of children */
    std::vector<WxHdf5NodeTypev*> children;
/** This node's type */
    int nodeType;
};

/**
 * WxHdf5IoTmpl does the reading and writing to HDF5 files for a given
 * data type.
 */
template <class DATATYPE>
class WxHdf5IoTmpl : public WxIoTmpl<DATATYPE> {

  public:

/**
 * Constructor.
 */
    WxHdf5IoTmpl();

/**
 * Virtual destructor
 */
    virtual ~WxHdf5IoTmpl();

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
        const std::vector<size_t>& dataSetLen, const DATATYPE* data) const;

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
        const std::vector<size_t>& dataSetLen, DATATYPE* data) const;

/**
 * Write an attribute.
 *
 * @param node the node to write under
 * @param attribName the name of the attribute
 * @param attrib the attribute to be written
 */
    virtual void writeAttribute(WxIoNodeType node,
        const std::string& attribName, const DATATYPE& attrib) const;

/**
 * Write an attribute.
 *
 * @param node the node to write under
 * @param attribName the name of the attribute
 * @param attrib the vector attribute to be written
 */
    virtual void writeVecAttribute(WxIoNodeType node,
        const std::string& attribName, const std::vector<DATATYPE>& attrib) const;

/**
 * Read an attribute.
 *
 * @param node the node to which this attribute belongs
 * @param attribName the name of the attribute
 * @param attrib the attribute to be read
 */
    virtual void readAttribute(WxIoNodeType node,
        const std::string& attribName, DATATYPE& attrib) const;

/**
 * Read an attribute.
 *
 * @param node the node to which this attribute belongs
 * @param attribName the name of the attribute
 * @param attrib the vector attribute to be read
 */
    virtual void readVecAttribute(WxIoNodeType node,
        const std::string& attribName, std::vector<DATATYPE>& attrib) const;

  private:

/** Private copy constructor to prevent use */
    WxHdf5IoTmpl(const WxHdf5IoTmpl<DATATYPE>&);

/** Private assignment to prevent use */
    WxHdf5IoTmpl<DATATYPE>& operator=(const WxHdf5IoTmpl<DATATYPE>&);

};

#endif // __wxhdf5iotmpl__
