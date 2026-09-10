"""Link against the conda environment the build is running inside.

SConstruct builds with a deliberately scrubbed environment:

    Environment(variables=vars, ENV={'PATH': ..., 'LM_LICENSE_FILE': ...,
                                     'LD_LIBRARY_PATH': ...})

so nothing the shell exported reaches the compiler or the linker. That is a
reasonable default - it stops a stray CFLAGS in someone's profile from changing
the build - but conda puts everything it needs in exactly those variables. Its
activation script sets

    LDFLAGS="-Wl,-rpath,$CONDA_PREFIX/lib -L$CONDA_PREFIX/lib ..."

and SCons discards it. The result is a build that compiles against the
environment's headers and then links against whatever the system happens to
have. On one cluster that produced a binary needing libgfortran.so.3 - a GCC 4
era library - from a conda GCC 15 toolchain whose own runtime is
libgfortran.so.5, so it linked and then would not start:

    error while loading shared libraries: libgfortran.so.3: cannot open
    shared object file: No such file or directory

Two things are needed and they are not the same. The -L decides which library
is chosen at LINK time; the rpath decides whether it can be found at RUN time.
Setting only the first gives a binary that builds and will not start unless
LD_LIBRARY_PATH happens to be right, and the SONAME recorded in it cannot be
fixed afterwards - it has to be relinked.

environment.yml in this repository tells people to build inside a conda
environment, so the build should behave as though it means it.
"""

import os


def conda_link_paths(environ=None):
    """(libpath, linkflags) additions for the active conda environment.

    Empty lists when not inside one, so this is inert on a distribution build.
    """
    environ = os.environ if environ is None else environ

    # CONDA_PREFIX is the active environment. PREFIX and BUILD_PREFIX are what
    # conda-build sets instead, where there is no activation at all.
    prefix = (environ.get('CONDA_PREFIX')
              or environ.get('PREFIX')
              or environ.get('BUILD_PREFIX'))
    if not prefix:
        return [], []

    lib = os.path.join(prefix, 'lib')
    if not os.path.isdir(lib):
        return [], []

    # -rpath rather than -rpath-link: the path has to survive into the binary,
    # because a batch job that runs the solver may not have activated anything.
    return [lib], ['-Wl,-rpath,' + lib]
