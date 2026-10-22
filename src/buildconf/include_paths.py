"""Keep the compiler's own system include directories out of -I.

`-I/usr/include` looks harmless and is not. The directory is already the last
entry in the compiler's search order, so naming it explicitly is redundant when
it is right - and when it is wrong it is fatal, because -I puts it AHEAD of the
toolchain's own sysroot headers.

That is how a conda build fails. conda-forge's GCC 15 makes _Float128, _Float32
and _Float64 built-in types and ships a matching glibc in its sysroot. An older
system glibc in /usr/include declares them as typedefs instead. Both are
self-consistent; the compiler only sees the conflict when -I/usr/include drags
the system copy in front of the sysroot:

    /usr/include/bits/floatn.h:86:20: error: redeclaration of C++ built-in
                                      type '_Float128' [-fpermissive]

Apollo put it there via this_host_config.py's `mpi_base = "/usr"`, which
config_mpi.py turns into -I$mpi_base/include. The same hazard applies to any
toolchain carrying its own libc headers: spack, a cross-compiler, or a module
built against a different sysroot.

Only the two directories the compiler already searches by default are dropped.
A path like /usr/include/hdf5/openmpi is a real location that is NOT on the
default path, so it stays.
"""

# The directories every hosted toolchain already searches. Anything below them
# is a genuine subdirectory and is left alone.
DEFAULT_SYSTEM_INCLUDES = ('/usr/include', '/usr/local/include')

# The same argument one stage later. -L/usr/lib puts the system libraries ahead
# of the toolchain's own, so a build that compiled against conda's headers can
# link the distribution's library instead - an ABI mismatch that shows up at
# run time rather than build time, which is worse than the compile error above.
# Only unambiguous defaults are listed: /usr/local/lib is deliberately absent
# because it is where a Homebrew or hand-built dependency legitimately lives.
DEFAULT_SYSTEM_LIBS = ('/usr/lib', '/usr/lib64', '/lib', '/lib64')


def _split(paths, defaults):
    kept, dropped = [], []
    for path in paths:
        normalised = str(path).rstrip('/') or '/'
        (dropped if normalised in defaults else kept).append(path)
    return kept, dropped


def drop_default_system_includes(paths):
    """(kept, dropped) - paths with the compiler's own defaults removed.

    Order is preserved. Trailing slashes are normalised so '/usr/include/' is
    recognised too.
    """
    return _split(paths, DEFAULT_SYSTEM_INCLUDES)


def drop_default_system_libs(paths):
    """(kept, dropped) - the same, for the linker's default directories."""
    return _split(paths, DEFAULT_SYSTEM_LIBS)
