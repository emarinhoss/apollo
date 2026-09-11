"""Where to look for Eigen3's headers.

Kept in a module of its own, and as a pure function of the environment, so it
can be tested. SConstruct is not importable, so anything decided inside it can
only be checked by running a build on a machine that happens to reproduce the
case - which is how this gap survived: the candidate list covered Debian,
Fedora, macOS and an explicit override, but not conda, while environment.yml in
this repository tells people to build inside a conda environment and installs
Eigen there.

Eigen ships its headers under a versioned directory. `#include <Eigen/Dense>`
therefore needs the *parent* of `Eigen` on the include path:

    $CONDA_PREFIX/include/eigen3/Eigen/Dense   ->  -I$CONDA_PREFIX/include/eigen3

Conda's activate script puts `$CONDA_PREFIX/include` on CPATH, which is why
PETSc and Boost are found without help and Eigen is not: they sit directly under
`include`, and Eigen sits one level deeper.
"""

import os
import sys


def eigen_candidates(environ=None, platform=None):
    """Directories to try, most specific first.

    `environ` and `platform` are injectable so the caller's real environment is
    not needed to test the decision.
    """
    environ = os.environ if environ is None else environ
    platform = sys.platform if platform is None else platform

    candidates = []

    # An explicit override wins over everything.
    override = environ.get('EIGEN3_INCLUDE_DIR')
    if override:
        candidates.append(override)

    # Conda, in the order a build is most likely to mean: the active
    # environment, then the prefixes conda-build sets.
    for key in ('CONDA_PREFIX', 'PREFIX', 'BUILD_PREFIX'):
        prefix = environ.get(key)
        if prefix:
            candidates.append(os.path.join(prefix, 'include', 'eigen3'))

    candidates += ['/usr/include/eigen3', '/usr/local/include/eigen3']

    if platform == 'darwin':
        for prefix in ('/opt/homebrew', '/usr/local'):
            candidates.append(os.path.join(prefix, 'include', 'eigen3'))

    # Preserve order, drop duplicates.
    seen = set()
    ordered = []
    for path in candidates:
        if path not in seen:
            seen.add(path)
            ordered.append(path)
    return ordered
