#include <apversion.h>

#include <petscversion.h>
#include <sstream>

namespace Apollo
{
  std::string
  version()
  {
    return APOLLO_GIT_DESCRIBE;
  }

  std::string
  buildInfo()
  {
    std::ostringstream os;
    os << "Apollo " << APOLLO_GIT_DESCRIBE
       << " (commit " << APOLLO_GIT_SHA << ")"
       << ", built " << __DATE__ << " " << __TIME__
       << ", PETSc " << PETSC_VERSION_MAJOR << "." << PETSC_VERSION_MINOR
       << "." << PETSC_VERSION_SUBMINOR;
#ifdef USE_BLAS
    os << ", BLAS";
#endif
#ifdef _OPENMP
    os << ", OpenMP";
#endif
#ifdef _DO_RANGE_CHECK_
    os << ", range-checked";
#endif
#if defined(__FAST_MATH__) || defined(APOLLO_FASTMATH)
    os << ", fast-math";
#endif
    // The single most important thing to know about a binary that produced a
    // suspicious result. Under -ffinite-math-only the compiler may assume no
    // value is ever NaN, which deletes the `x != x` checks the Euler and
    // nodal-DG kernels use to stop a diverged run - so a log without this
    // warning is a log whose silence about NaN means something.
#if defined(__FINITE_MATH_ONLY__) && __FINITE_MATH_ONLY__
    os << ", *** finite-math-only: NaN CHECKS DISABLED ***";
#endif
    return os.str();
  }
}
