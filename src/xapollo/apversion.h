#ifndef AP_VERSION_H
#define AP_VERSION_H

#include <string>

/**
 * Build provenance.
 *
 * A result produced by Apollo is only reproducible if you can tell which Apollo
 * produced it. APOLLO_GIT_DESCRIBE and APOLLO_GIT_SHA are defined on the compiler
 * command line by SConstruct, which reads them from git at configure time; they
 * fall back to "unknown" when the source is not a git checkout (a release
 * tarball, say).
 */

#ifndef APOLLO_GIT_DESCRIBE
#  define APOLLO_GIT_DESCRIBE "unknown"
#endif

#ifndef APOLLO_GIT_SHA
#  define APOLLO_GIT_SHA "unknown"
#endif

namespace Apollo
{
  /** Human-readable version, e.g. "v1.2-3-gabc1234" or "abc1234 (dirty)". */
  std::string version();

  /** One line naming the version, the build type, and the PETSc it links. */
  std::string buildInfo();
}

#endif // AP_VERSION_H
