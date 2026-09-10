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

Three files exist so you do not have to assemble this by hand:

| file | for |
| --- | --- |
| [`requirements.txt`](requirements.txt) | the Python side: SCons to build, numpy to analyse |
| [`requirements-dev.txt`](requirements-dev.txt) | the above plus linting and the legacy visualisation scripts |
| [`environment.yml`](environment.yml) | a conda environment that builds Apollo without root, for clusters with no PETSc module |

The C and C++ dependencies below still come from a package manager, a module
system, or that conda environment — `requirements.txt` cannot install PETSc.

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

That is the full set: on Ubuntu 22.04 and 24.04 these packages are enough to
build and run Apollo with no source builds and no edits to
`src/this_host_config.py`. It is also exactly what CI installs on both of those
releases, so it is checked on every push.

`libpetsc-real-dev` installs PETSc under `/usr/lib/petsc`, which the build finds
either automatically or via `export PETSC_DIR=/usr/lib/petsc`.

**Which PETSc you get depends on the release, and it matters.** Ubuntu 22.04
ships PETSc 3.15, 24.04 ships 3.19. Apollo builds against **3.11 and newer**:
`src/lib/petsc_compat.h` substitutes for the plex calls upstream has removed
over that range, and `test/test_petsc_compat.py` checks the substitutions fire
at the right versions — in both directions, so a shim cannot quietly collide
with a declaration that is still there. A PETSc older than 3.11 is untested and
will probably not compile.

If a build fails with a wall of template errors naming `DMPlex...` functions,
check the version first:

```bash
sed -n 's/#define PETSC_VERSION_\(MAJOR\|MINOR\|SUBMINOR\) *//p' \
    /usr/lib/petsc/include/petscversion.h | paste -sd. -
```

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

### Option 4: conda, for a cluster where you cannot install packages

```bash
conda env create -f environment.yml     # or mamba, which is much faster
conda activate apollo
cd src && scons build-opt
```

Read the comments in [`environment.yml`](environment.yml) first. If your cluster
provides PETSc, MPI and HDF5 as **modules**, prefer them — they are built against
the interconnect and conda-forge's OpenMPI is not — and use
`pip install -r requirements.txt` for the Python side only.

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

**A note on OpenMP.** It is off by default, for a correctness reason: the one
annotated loop (`subsolvers/wxnodaldg2dmethod.cc`) calls `DMPlexGetConeSize` and
`DMPlexPointLocalRead`/`Ref` once per element, and PETSc is not thread-safe
unless built `--with-threadsafety`, which the distribution packages are not.
Scaling is poor and non-monotonic besides — see
[Performance](#performance) for measured numbers. `openmp=yes` re-enables it.

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

> **Run one rank for anything you intend to analyse.** Earlier versions of this
> file gave the reason as incomplete output — one `.vtu` holding roughly 1/N of
> the cells. That was a bug in `test/vtu.py`, which flattened PETSc's per-rank
> `<Piece>` elements and kept only the last; the file was always complete, and a
> two-rank frame now reads back 4304 cells against the one-rank file's 4304.
> Fixing the reader exposed the actual problem: **the solver's answer depends on
> the rank count**. One rank and two give materially different fields on the
> Maxwell circular pulse as well as on the RMF decks, diverging from an identical
> initial frame. See [`docs/known-issues.md`](docs/known-issues.md) §15. The
> solver scales 3.88× on four ranks; that speed is not usable until this is
> understood.

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

**OpenMP is off by default**, for a correctness reason rather than a speed one.
Only one loop in the codebase carries an OpenMP directive, and it calls PETSc
(`DMPlexGetConeSize`, `DMPlexPointLocalRead`/`Ref`) once per element. PETSc is
not thread-safe unless it was configured `--with-threadsafety`, which
distribution packages are not - check for `PETSC_HAVE_THREADSAFETY` in your
`petscconf.h`.

Scaling is also poor and non-monotonic, because with PETSc's event logging on
(the default) each of those per-element calls updates shared counters. On the
bundled isentropic vortex, on a 4-core machine, best of three runs:

| `OMP_NUM_THREADS` | time | vs. serial |
| --- | --- | --- |
| 1 | 2.64 s | — |
| 2 | 1.68 s | 1.57x faster |
| 4 | 4.35 s | 1.65x slower |

Results are bit-identical across thread counts, so the parallelisation is
correct as far as Apollo's own state goes; it is PETSc's that is the problem.
Build with `openmp=yes` if you have a thread-safe PETSc and want to experiment,
and measure on your own hardware - four threads on four cores is a small sample.

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

All of these are run from the **repository root**, not from `src/`. Run them
from `src/` and `unittest discover` prints `Ran 0 tests` / `OK` and exits 0 —
a green run that tested nothing.

```bash
# The C++ unit tests: 86 checks across five binaries. Needs PETSc, not a
# built solver. Inside a conda environment, PETSC_DIR is $CONDA_PREFIX.
PETSC_DIR=/usr/lib/petsc make -C test/cxx

# Everything Python: tooling plus solver verification
python3 -m unittest discover -s test -v

# Just the parts that need no compiler
python3 -m unittest discover -s test -p 'test_python*' -v

# Solver against its own examples: needs a built binary
test/run_examples.sh                       # uses src/build-opt/apollo
test/run_examples.sh -b src/build-debug/apollo
test/run_examples.sh euler-isentropic-vortex
test/run_examples.sh -j 2                  # under mpirun - timing only, see below
```

**The solver-dependent tests skip themselves when no solver is present**, so a
green run does not by itself mean they ran. Check the skip count.

**`-j 2` is for timing, not for answers.** The solver's result depends on the
number of MPI ranks — one rank and two produce materially different fields
([`docs/known-issues.md`](docs/known-issues.md) §15, which carries a repro).
Use it to measure throughput; do not analyse what it writes.

There are two kinds of test here.

**Verification** (`test/test_vortex_accuracy.py`) asks whether Apollo solves the
equations it claims to. The isentropic vortex is an exact, smooth solution of the
Euler equations that translates without changing shape, so the difference from
the analytic result at t=1 is pure discretization error. The suite measures it
for every numerical flux, and measures the convergence rate across the 290/1110/
4454-cell meshes shipped with the example:

```
    grid.msh       290 cells   L2(rho) = 1.58e-02
    grid2.msh     1110 cells   L2(rho) = 5.13e-03      observed order 1.68
    grid3.msh     4454 cells   L2(rho) = 1.14e-03      observed order 2.16
```

Second order is what P1 DG should give. This is the test that catches a wrong
answer, as opposed to a crash: a flux that writes the wrong momentum component
runs to completion and returns nonsense, and only a comparison against a known
solution notices.

**Smoke** (`test/run_examples.sh` and the rest of `test/`) asks whether things
still run at all.

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

## Known issues

[docs/known-issues.md](docs/known-issues.md) records defects that are confirmed
but not fixed, with the evidence and a way to reproduce each. Read it before
trusting a result from the multifluid module, before using
`Numerical_Flux = Wave`, before enabling a slope limiter, or before assuming a
run can be restarted.

[docs/rmf-frc-model-assessment.md](docs/rmf-frc-model-assessment.md) assesses the
`multifluid/rmf_frc` example against the RMF current-drive literature: what is
right (the model class, the regime, the formation scenario, the flux conserver),
what is not, and a phased plan with a verification for each step. Three phases
of that plan are done — the friction energetics; moving the RMF drive from the
plasma edge to an antenna current so that the plasma can screen and load it
(`examples/unstructuredDG/multifluid/rmf_frc/antenna/`, verified against a closed
form by `test/test_rmf_antenna_field.py`); and characteristic-consistent
injection at the edge-driven boundary, which imposes only the incoming Riemann
invariants and takes the outgoing ones from the interior
(`src/hyperapps/maxwell/apmaxwellcharacteristics.h`, opt-in per deck via `c0`).
Phase 3, the validation campaign against the literature, is partly done: the
question of which way each drive turns is settled (`test/cxx/`, 15 checks), and
the diagnostics and scan harness the remaining items need are built and verified
against closed forms (`scripts/rmf_diagnostics.py`, `scripts/rmf_scan.py`). The
comparisons themselves are not run: the assessment says what each would cost —
measured, not guessed — and records two findings that block them, both about the
shipped decks rather than the code. The reduced speed of light is 1.2% above the
electron sound speed and below the fastest electron characteristic of the
penetrated state the campaign exists to measure; raising it trades away Debye
resolution one for one, since the product of the two is fixed by the density, and
the deck is already under-resolved there. And the edge-driven boundary condition
drives the divergence-cleaning potential without bound, so within a few hundred
steps its transverse field depends on a numerical parameter — the antenna deck
does not do this. The assessment also records why the antenna deck has no vacuum
region around its plasma column, which is a limit of the model rather than a
choice.

The runs Phase 3 asks for are laid out ready to start in
[`examples/unstructuredDG/multifluid/rmf_frc/phase3/`](examples/unstructuredDG/multifluid/rmf_frc/phase3/),
one self-contained folder each, with a README saying what to run, in what order,
and what it costs. Start with `00-divergence-gate`: it takes six minutes and
decides whether the other three folders are worth their compute. There is a
SLURM job-array template and a Colab bootstrap script beside them.

Four tools support that work, and each is verified against closed forms by tests
that need no solver:

| script | what it does |
| --- | --- |
| `scripts/rmf_diagnostics.py` | what a run produced: driven field, ζ, penetration, current-layer thickness |
| `scripts/rmf_scan.py` | parameter scans, and what one would cost before you start it |
| `scripts/rmf_cleaning_check.py` | the gate: does the answer depend on the divergence-cleaning speed? |
| `scripts/rmf_make_phase3_runs.py` | regenerates the run folders from the shipped decks |

## Documentation

Documentation is a work in progress. The most reliable sources are, in order:
- this file, for building and running;
- [CLAUDE.md](CLAUDE.md), for the conventions and the traps worth knowing before
  editing anything;
- `examples/unstructuredDG/`, for what a deck looks like per physics module;
- `test/run_examples.sh`, for a working end-to-end invocation;
- [docs/known-issues.md](docs/known-issues.md), for what is broken;
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