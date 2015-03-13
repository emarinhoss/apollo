#ifndef WXUBLASMATVEC_H
#define WXUBLASMATVEC_H

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

namespace numeric_containers
{
    namespace ublas = boost::numeric::ublas; // basically boost ublas will be the implementation for basic matrix container

    typedef ublas::matrix<double> DMAT; // matrix of doubles
    typedef ublas::matrix<int>    IMAT; // matrix of int

    typedef ublas::vector<double> DVec; // vector of doubles
    typedef ublas::vector<int> IVec;    // vector of int

}

#endif // WXUBLASMATVEC_H
