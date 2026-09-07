# Apollo

A high-performance Discontinuous Galerkin (DG) finite element solver for hyperbolic partial differential equations (PDEs) on unstructured grids. Apollo is designed for plasma physics simulations, nuclear fusion research, and computational fluid dynamics applications.

## Overview

Apollo implements a flexible and extensible DG framework capable of solving multiple physics systems including electromagnetic wave propagation, compressible fluid dynamics, and multi-species plasma transport. The solver leverages PETSc for massively parallel computations and supports high-order spatial discretizations on complex geometries.

## Features

- **Multiple Physics Modules:**
  - Maxwell equations (electromagnetics, RF wave propagation)
  - Euler equations (compressible gas dynamics)
  - Multifluid models (plasma physics with multiple species)
  - Advection equations

- **Advanced Numerics:**
  - High-order Discontinuous Galerkin method
  - Unstructured mesh support (2D/3D)
  - Entropy-based formulations for stability
  - Shock capturing and limiters

- **Scalable Parallelization:**
  - MPI-based domain decomposition via PETSc
  - Efficient sparse matrix operations
  - Scalable to large-scale HPC systems

- **Application Areas:**
  - Nuclear fusion devices (FRC, RMF heating systems)
  - Space propulsion and plasma thrusters
  - Computational electromagnetics
  - Hypersonic flow simulations

## Requirements

### Core Dependencies

- **C++ compiler** - GCC or Clang with C++14 support (GCC 5+, Clang 3.4+).
  The build passes `-std=c++14`.
- **MPI library** - OpenMPI 1.8+ or MPICH 3.0+. The build uses the `mpicc` and
  `mpicxx` wrappers.
- **PETSc** - 3.6 or later, built with MPI. Regularly built against 3.19.
  Apollo still calls a few routines PETSc deprecated in 3.8 (`TSSetDuration`,
  `TSSetInitialTimeStep`, `TSGetTimeStepNumber`), which compile with warnings.
- **SCons** - the build system (`pip install scons`). Apollo does **not** use
  CMake.
- **Python 3** - required at build time (SCons) and to preprocess input decks.

### Additional Libraries

- **HDF5** - required in practice. Apollo does not call HDF5 itself, but PETSc's
  `petscviewerhdf5.h` includes `<hdf5.h>`, so any HDF5-enabled PETSc (which
  includes every distribution package) makes its headers a build requirement.
- **Boost** - required by the multifluid module (`format`, uBLAS, Graph).
  Headers only; no Boost libraries are linked.
- **Eigen 3** - required by the multifluid module (`<Eigen/Dense>`).
- **GSL** (GNU Scientific Library) - required; `libgsl` and `libgslcblas` are
  linked unconditionally.
- **BLAS** - optional but recommended. When a BLAS with `<cblas.h>` is found,
  the cubature matrix products go through `cblas_dgemm`; otherwise a built-in
  loop is used. Both build variants follow the same choice.
- **Exodus II** - optional, for Exodus mesh I/O. The build warns and continues
  without it.

## Installation

### Option 1: Install Dependencies via System Package Manager

#### Ubuntu/Debian

```bash
sudo apt-get update
sudo apt-get install build-essential git wget python3-pip
sudo apt-get install libopenmpi-dev openmpi-bin
sudo apt-get install libpetsc-real-dev      # PETSc, no source build needed
sudo apt-get install libhdf5-openmpi-dev
sudo apt-get install libboost-dev
sudo apt-get install libeigen3-dev
sudo apt-get install libgsl-dev
sudo apt-get install libopenblas-dev        # optimized BLAS, provides <cblas.h>
pip3 install scons
```

That is the full set: on Ubuntu 24.04 these packages are enough to build and run
Apollo with no source builds and no edits to `src/this_host_config.py`. It is
also exactly what CI installs, so it is checked on every push.

`libpetsc-real-dev` installs PETSc under `/usr/lib/petsc`, which the build finds
either automatically or via `export PETSC_DIR=/usr/lib/petsc`.

#### CentOS/RHEL/Fedora

```bash
sudo yum install gcc gcc-c++ git wget python3-pip
sudo yum install openmpi openmpi-devel
sudo yum install hdf5-openmpi hdf5-openmpi-devel
sudo yum install boost boost-devel
sudo yum install gsl gsl-devel
sudo yum install openblas openblas-devel  # Optimized BLAS for better performance
sudo yum install blas blas-devel lapack lapack-devel
```

Load MPI module:
```bash
module load mpi/openmpi-x86_64
```

#### macOS

Using [Homebrew](https://brew.sh/) package manager:

```bash
# Install Homebrew if not already installed
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# Install dependencies
brew install gcc git wget
brew install open-mpi
brew install hdf5-mpi
brew install boost
brew install gsl
brew install openblas  # Optimized BLAS for better performance
brew install lapack
brew install libomp    # OpenMP support for Clang (required for parallel performance)
brew install eigen     # C++ template library for linear algebra
brew install seacas    # For Exodus II mesh I/O (optional but recommended)
```

**Note for macOS users:**
- Homebrew installs libraries in `/usr/local` (Intel) or `/opt/homebrew` (Apple Silicon)
- You may need to specify compiler paths when building Apollo:
  ```bash
  export CC=/usr/local/bin/gcc-13
  export CXX=/usr/local/bin/g++-13
  ```
- For Apple Silicon Macs, add Homebrew to your PATH:
  ```bash
  export PATH=/opt/homebrew/bin:$PATH
  ```

For PETSc on macOS, you can either:
- Install via Homebrew: `brew install petsc` (easiest option)
- Build from source using Option 3 below (for custom configuration)

**Exodus II library:**
- The `seacas` package provides the Exodus II mesh I/O library
- If `seacas` is not available via Homebrew, you can build PETSc with built-in Exodus support (see Option 3 with `--download-exodusii=1`)
- Apollo can build without Exodus but mesh I/O capabilities will be limited

### Option 2: Build Dependencies from Source

Prefer Option 1 where the packages exist; only PETSc is likely to need building,
and only if your distribution has no suitable package.

The `makefile` in the repository root is an unfinished dependency bootstrapper
and does not work: `make all` names three targets that do not exist
(`builddep_hdf5`, `builddep_boost`, `builddep_uc`), and its one real target
expands `MPINAME`, `MPIZIP` and `MPIINSTALLDIR`, none of which are defined. Do
not use it. It is kept only as a starting point for a future bootstrap script.

### Option 3: Install PETSc Manually

PETSc is a critical dependency and often needs to be built from source:

```bash
# Download PETSc
cd $HOME/software
wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.6.0.tar.gz
tar -xzf petsc-3.6.0.tar.gz
cd petsc-3.6.0

# Configure PETSc with required options
./configure --prefix=$HOME/software/petsc \
  --with-mpi=1 \
  --download-fblaslapack=1 \
  --download-exodusii=1 \
  --with-hdf5=1 \
  --with-debugging=0 \
  COPTFLAGS='-O3' \
  CXXOPTFLAGS='-O3' \
  FOPTFLAGS='-O3'

# Build and install
make PETSC_DIR=$HOME/software/petsc-3.6.0 PETSC_ARCH=arch-linux-c-opt all
make PETSC_DIR=$HOME/software/petsc-3.6.0 PETSC_ARCH=arch-linux-c-opt install
```

Set environment variables:
```bash
export PETSC_DIR=$HOME/software/petsc
export PETSC_ARCH=arch-linux-c-opt
```

### Building Apollo

Apollo uses SCons as its build system. Once all dependencies are installed:

```bash
cd /path/to/apollo/src

# Build optimized version (recommended for production)
scons build-opt

# Or build debug version (for development)
scons build-debug
```

This produces `src/build-opt/apollo` (or `src/build-debug/apollo`). Note the
name: the executable is `apollo`.

**If a dependency is not found**, the build stops and tells you what to install
and how to point it at an existing installation. The search order for PETSc is:
`petsc_base` on the command line, then `src/this_host_config.py`, then
`$PETSC_DIR` (honouring `$PETSC_ARCH`), then the compiler's default paths.

**Build variables.** Anything set in `src/this_host_config.py` can be overridden
on the command line, which always wins. Run `scons -h` for the full list.

| Variable | Default | Purpose |
| --- | --- | --- |
| `petsc_base` | *(autodetect)* | PETSc installation prefix |
| `mpi_base` | `/usr` | MPI installation prefix |
| `hdf5_base` | *(autodetect)* | HDF5 installation prefix |
| `boost_base`, `blas_base` | *(system paths)* | Boost / BLAS prefixes |
| `arch` | `native` | value for `-march=`; use a baseline such as `x86-64-v2`, or `none`, to build a portable binary |
| `tune` | *(follows `arch`)* | value for `-mtune=`; only CPU names, not ISA levels |
| `fastmath` | `no` | add `-ffast-math` (see below) |
| `openmp` | `no` | add `-fopenmp` (see below) |

```bash
# Build for a cluster whose compute nodes are older than the login node
scons build-opt arch=x86-64-v2

# Point at a PETSc that is not on the default path
scons build-opt petsc_base=$HOME/software/petsc
```

**A note on `-ffast-math`.** It is not enabled by default, and turning it on is
not free. `-ffast-math` implies `-ffinite-math-only`, which permits the compiler
to assume no value is ever NaN and therefore to fold `x != x` to `false` - and
`x != x` is the idiom Apollo uses to detect a diverged solution in the Euler
equations and the nodal-DG element loop. With `-ffast-math` those guards are
deleted from the binary and a run that has gone to NaN keeps going silently.
`fastmath=yes` adds `-fno-finite-math-only` alongside it so the guards survive,
but the remaining reassociation still changes results run to run.

**A note on OpenMP.** It is off by default. Only one loop is annotated
(`subsolvers/wxnodaldg2dmethod.cc`), and that loop calls `DMPlexGetConeSize` and
`DMPlexPointLocalRead`/`Ref` once per element. PETSc is not thread-safe unless
built `--with-threadsafety`, which the distribution packages are not - check for
`PETSC_HAVE_THREADSAFETY` in your `petscconf.h`. With PETSc's logging enabled
(the default) every such call also updates shared counters, and threaded runs
measured slower than serial rather than faster. `openmp=yes` re-enables it if
you have a thread-safe PETSc and want to experiment.

**What the optimized build actually does:**

- `-O3 -funroll-loops -ftree-vectorize -fno-math-errno`, plus `-march`/`-mtune`
  per the `arch` and `tune` variables (`native` by default).
- `cblas_dgemm` for the cubature matrix products when a BLAS with `<cblas.h>` is
  present. `build-debug` makes the same choice, so the two variants agree
  numerically.
- No `-ffast-math` and no OpenMP unless you ask for them; see the notes above.

To get the most out of it: install OpenBLAS (`sudo apt-get install
libopenblas-dev`) and scale with MPI ranks (`mpirun -np N`), which is where
Apollo's parallelism actually lives.

**Note:** If dependencies are installed in non-standard locations, you may need to set environment variables:

```bash
export PATH=$HOME/apollo/software/bin:$PATH
export LD_LIBRARY_PATH=$HOME/apollo/software/lib:$LD_LIBRARY_PATH
export PETSC_DIR=$HOME/apollo/software/petsc
```

## Usage

### Input decks: `.pin` and `.inp`

Running a simulation is two steps, because the deck you edit is not the deck the
solver reads.

The `.pin` files under `examples/` are **templates**: Python that runs at
preprocessing time to compute derived quantities, so a deck can say
`cfl = 1.0/3.0` or `PI = math.pi` instead of hard-coding a number.
`scripts/wxinpparse.py` macro-expands and preprocesses a `.pin` into a plain
`.inp`, and `.inp` is what the solver parses. `.inp` files are generated and are
therefore listed in `.gitignore`.

```bash
# 1. Expand the deck (writes isentropicVortex.inp next to the .pin)
PYTHONPATH=scripts python3 scripts/wxinpparse.py -i isentropicVortex.pin

# 2. Run
src/build-opt/apollo -i isentropicVortex.inp
```

Handing the solver a `.pin` directly gives you a parse error on the first Python
line of the file.

### Basic execution

The executable is **`apollo`**, and it takes the input file with **`-i`**.
A bare filename argument is ignored, and Apollo then looks for its default,
`apollo.inp`.

```bash
# Serial
src/build-opt/apollo -i <input>.inp

# Parallel
mpirun -np <ranks> src/build-opt/apollo -i <input>.inp
```

Mesh files named by a deck are resolved relative to the **working directory**,
not to the deck, so run from the directory holding the mesh (or copy the mesh in).

Options, as reported by `apollo --help`:

| Option | Meaning |
| --- | --- |
| `-i <file>`, `--input-file=<file>` | input file; defaults to `apollo.inp` |
| `-o <prefix>`, `--output-prefix=<prefix>` | output prefix; defaults to the input name without its extension |
| `-r <file>`, `--restart=<file>` | restart from `<file>` |
| `--help` | print the option summary |

Two entries in that summary do not do what they say: the short form `-h` is not
accepted (use `--help`), and `--real-type=float` is parsed but has no effect,
since `main()` instantiates the solver as `double` and the `float`
instantiations are commented out.

### Example simulations

Each example directory holds a `.pin` deck and the `.msh` meshes it references.

```bash
cd examples/unstructuredDG/euler/isentropicVortex
PYTHONPATH=../../../../scripts python3 ../../../../scripts/wxinpparse.py -i isentropicVortex.pin
../../../../src/build-opt/apollo -i isentropicVortex.inp
```

The same pattern applies to `unstructuredDG/advection/advection.pin`,
`unstructuredDG/maxwell/circularPulse/circularPulse.pin` and the multifluid
decks under `unstructuredDG/multifluid/`.

`test/run_examples.sh` does all of this for the three quickest cases; see
[Testing](#testing).

## Example Problems

The repository includes various benchmark and application problems:

- **Maxwell:** Circular pulse, transverse magnetic modes, RMF antenna modeling for FRC
- **Euler:** Isentropic vortex, forward/backward facing steps, scramjet inlet, explosion tests
- **Multifluid:** Collisional-radiative models, RMF-driven FRC plasmas

Visualization scripts (Python) are provided in the `scripts/` directory for post-processing results.

## Performance

Apollo's parallelism is MPI domain decomposition through PETSc. That is the axis
that scales; use `mpirun -np N`.

### What the optimized build gives you

- **Compiler optimization** - `-O3`, loop unrolling, auto-vectorization, and
  `-march`/`-mtune` (`native` by default, so the build machine's SIMD width is
  used). Set `arch=` to a baseline such as `x86-64-v2` when the binary has to run
  on nodes older than the one you built on, or it will die with SIGILL.
- **BLAS** - `cblas_dgemm` replaces the hand-written triple loop for the cubature
  matrix products when a BLAS with `<cblas.h>` is available. Install OpenBLAS.

### What it does not give you

**OpenMP is off by default and should stay off** unless you have specifically
checked otherwise. Only one loop in the codebase carries an OpenMP directive, and
it calls PETSc (`DMPlexGetConeSize`, `DMPlexPointLocalRead`/`Ref`) once per
element. PETSc is not thread-safe unless it was configured
`--with-threadsafety`, which distribution packages are not; and with PETSc's
event logging on - the default - each of those calls updates shared counters, so
adding threads adds contention. Threaded runs of the bundled examples measured
slower than serial. See the `openmp` build variable if you want to experiment.

**`-ffast-math` is off by default.** See the note under [Building
Apollo](#building-apollo): it deletes Apollo's own NaN guards.

Earlier revisions of this file advertised "10-60x" from OpenMP plus BLAS plus
fast-math. That figure was never measured and the OpenMP part of it was
negative; it has been removed rather than restated.

### Measuring

```bash
# PETSc's own profiler: where the time goes, per stage and per event
PETSC_OPTIONS="-log_view" src/build-opt/apollo -i input.inp

perf stat -e cycles,instructions,cache-misses src/build-opt/apollo -i input.inp
```

### Known optimization opportunities

Unimplemented, listed roughly by expected value:

- Overlap halo exchange with interior element work.
- Make the threaded element loop viable: hoist the per-element PETSc lookups out
  of the loop so the parallel region touches no PETSc state, then re-enable
  OpenMP.
- Reorder elements for locality.
- GPU offload of the element kernels.

## Testing

```bash
# Python tooling: ~1 s, needs only python3 (pyflakes adds one more check)
python3 -m unittest discover -s test -v

# Solver against its own examples: ~15 s, needs a built binary
test/run_examples.sh                       # uses src/build-opt/apollo
test/run_examples.sh -b src/build-debug/apollo
test/run_examples.sh -j 2                  # under mpirun
test/run_examples.sh euler-isentropic-vortex
```

Or through the root makefile: `make test`, `make test-python`,
`make test-examples`; `make help` lists everything, including `make opt`,
`make debug` and `make deps-ubuntu`.

`test/run_examples.sh` preprocesses each deck and runs the solver to completion,
failing on a non-zero exit or on a NaN, negative-pressure or exception
diagnostic in the log - the solver does not always exit non-zero on those.

Both suites, plus a syntax check over every SCons script and a full build of
both variants against distribution packages, run in CI on every push
(`.github/workflows/ci.yml`).

## About PETSc

Apollo uses PETSc (Portable, Extensible Toolkit for Scientific Computation) for parallel numerical computations. PETSc provides:

- Highly portable parallel framework (distributed memory, shared memory, hybrid MPI-GPU)
- Comprehensive suite of parallel linear/nonlinear solvers
- Scalability to millions/billions of unknowns
- Industry-standard MPI-based communication
- Extensive documentation and active development community

## Reproducibility

Every run begins by recording which binary produced it:

```
Apollo v1.0-3-g5af1384 (commit 5af13846e732), built Sep  7 2026 19:14:18, PETSc 3.19.6, BLAS
Input file: isentropicVortex.inp
MPI ranks: 4
```

The version comes from `git describe` at build time, with `-dirty` appended when
the working tree has uncommitted changes, so a log can be traced back to a
commit. The same line is printed by `apollo --help`.

The tail of that line names the build options that change results: `BLAS`,
`OpenMP`, `range-checked`, `fast-math`, and

```
*** finite-math-only: NaN CHECKS DISABLED ***
```

which appears when the binary was compiled such that its NaN guards were
optimized away. If you are looking at a suspicious result, check that line first.

## Repository layout

| Path | Contents |
| --- | --- |
| `src/SConstruct`, `src/*/SConscript` | build definition |
| `src/buildconf/` | dependency detection (MPI, PETSc, HDF5, BLAS, Boost) |
| `src/this_host_config.py` | host-specific build defaults; command line overrides it |
| `src/xapollo/` | `main()` and the simulation driver |
| `src/lib/` | array and matrix types, nodal-DG operators, quadrature tables, input-deck parser, logging, MPI wrappers |
| `src/hyper/` | hyperbolic equation and source-term base classes |
| `src/hyperapps/` | physics: `advection/`, `euler/`, `maxwell/`, `multifluid/` |
| `src/solvers/`, `src/subsolvers/` | time stepping, limiters, boundary conditions, the DG element loop |
| `src/lapack_lite/` | vendored LAPACK subset, **not currently built** |
| `scripts/` | input preprocessing and post-processing / visualization |
| `examples/` | `.pin` decks and `.msh` meshes per physics module |
| `test/` | Python regression tests and the solver smoke tests |

Two naming conventions coexist: `wx`-prefixed files are inherited from WarpM, the
project Apollo descends from, and `ap`-prefixed files are Apollo's own. They are
the same codebase; the prefix carries no meaning beyond age.

## Post-processing

`scripts/` holds the visualization tooling. Beyond `numpy`, individual scripts
import `matplotlib`/`pylab`, `pyvtk`, `tvtk`/`mayavi`, `tables` (PyTables) and
`joblib`; install what the script you need asks for. `wxdata.py` is the common
reader, and `wxplot.py` a command-line front end to it. The `wxxdmf*.py` family
writes XDMF wrappers for the HDF5 output so it opens in ParaView or VisIt.

Be aware that several scripts are near-duplicates of each other
(`wxdata.py`/`wxdata_3949.py`, `wxunsdgdata.py`/`wxunsdgdata2.py`,
`SGC.py`/`SGC_MB.py`, and the `wxxdmf*` set), and it is not always obvious which
is canonical.

## Documentation

Documentation is a work in progress. The most reliable sources are, in order:
- this file, for building and running;
- `examples/unstructuredDG/`, for what a deck looks like per physics module;
- `test/run_examples.sh`, for a working end-to-end invocation;
- the headers in `src/`, for the class structure.

There is no generated API documentation and no description of the discretization
beyond the source.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md).

## Citing Apollo

See [CITATION.cff](CITATION.cff), or use GitHub's "Cite this repository" button.

## Author

Created by Eder Sousa