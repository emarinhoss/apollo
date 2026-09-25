// std includes

// hdf5 includes
#define H5_USE_16_API
#include <hdf5.h>

// WarpX lib includes
#include "wxhdf5io.h"
#include "wxhdf5iotmpl.h"
#include "wxexcept.h"

WxHdf5Io::WxHdf5Io(MPI_Comm mc, MPI_Info mi) : WxIoBase(".h5") {
// Store mpi stuff, create templated IO objects
  setup(mc, mi);
}

WxHdf5Io::WxHdf5Io(const std::string& bn, int d, MPI_Comm mc,
  MPI_Info mi) : WxIoBase(bn, d, ".h5") {
// Store mpi stuff, create templated IO objects
  setup(mc, mi);
}

WxHdf5Io::~WxHdf5Io() {
// Must close any open files.
  closeOpenFiles();
}


void WxHdf5Io::setup(MPI_Comm mc, MPI_Info mi) {

// Store mpi stuff
  mpiComm = mc;
  mpiInfo = mi;

// Create the templated io objects
  this->addIo( new WxHdf5IoTmpl<char>());
  this->addIo( new WxHdf5IoTmpl<unsigned char>() );
  this->addIo( new WxHdf5IoTmpl<short>() );
  this->addIo( new WxHdf5IoTmpl<unsigned short>() );
  this->addIo( new WxHdf5IoTmpl<int>() );
  this->addIo( new WxHdf5IoTmpl<unsigned int>() );
  this->addIo( new WxHdf5IoTmpl<long>() );
  this->addIo( new WxHdf5IoTmpl<unsigned long>() );
  this->addIo( new WxHdf5IoTmpl<long long>() );
  //this->addIo( new WxHdf5IoTmpl<float>() );
  this->addIo( new WxHdf5IoTmpl<double>() );
  this->addIo( new WxHdf5IoTmpl<long double>() );
}

WxIoNodeType WxHdf5Io::createFile(const std::string& fileName) {
// Determine access properties
  hid_t plistId = H5Pcreate(H5P_FILE_ACCESS);
#ifdef _DO_USE_MPI_
  H5Pset_fapl_mpio(plistId, mpiComm, mpiInfo);
#endif
  long fn = H5Fcreate(fileName.c_str(), H5F_ACC_TRUNC,
    H5P_DEFAULT, plistId);
  WxIoNodeTypev *fileNode = new WxHdf5NodeTypev(fn);
  H5Pclose(plistId);
  if (fn < 0) {
    WxExcept wxe("WxHdf5Io::WxHdf5Io: unable to create file ");
    wxe << "'" << fileName << "'";
    throw wxe;
  }
  addOpenFile(fileNode);
  return fileNode;
}

WxIoNodeType WxHdf5Io::openFile(const std::string& fileName,
  const std::string& perms) {
// Determine MPI access properties
  hid_t plistId = H5Pcreate(H5P_FILE_ACCESS);
#ifdef _DO_USE_MPI_
  H5Pset_fapl_mpio(plistId, mpiComm, mpiInfo);
#endif
// Determine whether writable
  long fn;

  if (perms.find_first_of('w') != std::string::npos) {
    fn = H5Fopen(fileName.c_str(), H5F_ACC_TRUNC, plistId);
  }
  else {
    fn = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, plistId);
  }
  WxIoNodeTypev *fileNode = new WxHdf5NodeTypev(fn);
  H5Pclose(plistId);
  if (fileNode < 0) {
    WxExcept wxe("WxHdf5Io::WxHdf5Io: unable to create file ");
    wxe << "'" << fileName << "'";
    throw wxe;
  }
  addOpenFile(fileNode);
  return fileNode;
}

WxIoNodeType
WxHdf5Io::createGroup(WxIoNodeType node, const std::string& dataName) const {
  long dn = H5Gcreate(
    static_cast<WxHdf5NodeTypev*>(node)->node,
    dataName.c_str(),
    0);
  WxIoNodeTypev *dataNode = new WxHdf5NodeTypev(dn, static_cast<WxHdf5NodeTypev*>(node));
  static_cast<WxHdf5NodeTypev*>(dataNode)->setNodeType(0);
  return dataNode;
}

WxIoNodeType
WxHdf5Io::openGroup(WxIoNodeType node, const std::string& dataName) const {
  long dn = H5Gopen(
    static_cast<WxHdf5NodeTypev*>(node)->node,
    dataName.c_str());
  WxIoNodeTypev* dataNode = new WxHdf5NodeTypev(dn, static_cast<WxHdf5NodeTypev*>(node));
  static_cast<WxHdf5NodeTypev*>(dataNode)->setNodeType(0);
  return dataNode;
}

WxIoNodeType
WxHdf5Io::createDataSet(WxIoNodeType node, const std::string& dataName) const {
  hid_t fileSpace = H5Screate(H5S_SCALAR);
  long dn = H5Dcreate(
    static_cast<WxHdf5NodeTypev*>(node)->node,
    dataName.c_str(),
    H5T_NATIVE_INT,
    fileSpace,
    H5P_DEFAULT);
  WxIoNodeTypev *dataNode = new WxHdf5NodeTypev(dn, static_cast<WxHdf5NodeTypev*>(node));
  H5Sclose(fileSpace);

  return dataNode;
}

WxIoNodeType
WxHdf5Io::openDataSet(WxIoNodeType node, const std::string& dataName) const {
  long dn = H5Dopen(
    static_cast<WxHdf5NodeTypev*>(node)->node,
    dataName.c_str());
  WxIoNodeTypev* dataNode = new WxHdf5NodeTypev(dn, static_cast<WxHdf5NodeTypev*>(node));

  return dataNode;
}

void WxHdf5Io::closeFile(WxIoNodeTypev *fileNode) {
  removeOpenFile(fileNode);	// Remove reference in base class
  delete fileNode;
}

void WxHdf5Io::closeDataSet(WxIoNodeTypev *node) const {
// Close a data set
  delete node;
}

void
WxHdf5Io::writeStrAttribute(WxIoNodeType node, const std::string& attribName, const std::string& attrib) const 
{
  // Add attribute that saves a vector of chars
  hsize_t size = attrib.size();
  hid_t attrDS = H5Screate(H5S_SCALAR);
  // Create attribute to store the number of elements
  hid_t attrTP = H5Tcopy(H5T_C_S1);
  H5Tset_size(attrTP, size);
  hid_t attrID = H5Acreate(static_cast<WxHdf5NodeTypev*>(node)->node,
    attribName.c_str(), attrTP, attrDS, H5P_DEFAULT);
  // Write out the attribute
  H5Awrite(attrID, attrTP, attrib.c_str());
  //this->comm()->barrier();
  // Close the attribute
  H5Aclose(attrID);
  // Close the attribute type
  H5Tclose(attrTP);
  // Close the attribute dataspace
  H5Sclose(attrDS);
}


