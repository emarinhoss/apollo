#ifndef __wxmpitraits__
#define __wxmpitraits__

#include <mpi.h>
#include "wxmsgbase.h"

/**
 * Type traits for use in MPI messengers
 */
template<typename T> class WxMpiTraits;

// char
template<> class WxMpiTraits<char> { 
  public:
    static MPI_Datatype mpiType() {
      return MPI_CHAR;
    };
};
// unsigned char
template<> class WxMpiTraits<unsigned char> { 
  public:
    static MPI_Datatype mpiType() {
      return MPI_UNSIGNED_CHAR;
    }; 
};
// short
template<> class WxMpiTraits<short> {
  public:
    static MPI_Datatype mpiType() {
      return MPI_SHORT;
    }; 
};
// unsigned short
template<> class WxMpiTraits<unsigned short> {
  public:
    static MPI_Datatype mpiType() {
      return MPI_UNSIGNED_SHORT;
    }; 
};
// int
template<> class WxMpiTraits<int> {
  public:
    static MPI_Datatype mpiType() {
      return MPI_INT;
    };
};
// unsigned
template<> class WxMpiTraits<unsigned> {
  public:
    static MPI_Datatype mpiType() {
      return MPI_UNSIGNED;
    };
};
// long
template<> class WxMpiTraits<long> {
  public:
    static MPI_Datatype mpiType() {
      return MPI_LONG;
    };
};
// unsigned long
template<> class WxMpiTraits<unsigned long> {
  public:
    static MPI_Datatype mpiType() {
      return MPI_UNSIGNED_LONG;
    };
};
// float
template<> class WxMpiTraits<float> {
  public:
    static MPI_Datatype mpiType() {
      return MPI_FLOAT;
    };
};
// double
template<> class WxMpiTraits<double> {
  public:
    static MPI_Datatype mpiType() {
      return MPI_DOUBLE;
    };
};
// long double
template<> class WxMpiTraits<long double> {
  public:
    static MPI_Datatype mpiType() {
      return MPI_LONG_DOUBLE;
    };
};
// long long int
template<> class WxMpiTraits<long long int> {
  public:
    static MPI_Datatype mpiType() {
      return MPI_LONG_LONG_INT;
    };
};

/**
 * Type traits for use in MPI all-reduce operators
 */
template<const unsigned T> class WxMpiAllReduce;

// MIN
template<> class WxMpiAllReduce<WX_MSG_MIN> { 
  public:
    static MPI_Op mpiOp() {
      return MPI_MIN;
    };
};
// MAX
template<> class WxMpiAllReduce<WX_MSG_MAX> { 
  public:
    static MPI_Op mpiOp() {
      return MPI_MAX;
    };
};
// SUM
template<> class WxMpiAllReduce<WX_MSG_SUM> {
  public:
    static MPI_Op mpiOp() {
      return MPI_SUM;
    };
};

#endif // __wxmpitraits__
