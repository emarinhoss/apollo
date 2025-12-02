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

- **C++ compiler** - GCC 4.8+ or Clang with C++11 support
- **MPI library** - OpenMPI 1.8+ or MPICH 3.0+
- **PETSc** - Version 3.6.0 or later (with MPI support)
- **CMake** - Version 2.8 or later

### Additional Libraries

- **HDF5** - For parallel I/O and data storage
- **Exodus II** - For unstructured mesh I/O
- **GSL** (GNU Scientific Library) - For numerical routines
- **Boost** - C++ utility libraries
- **BLAS/LAPACK** - For linear algebra operations

## Installation

### Option 1: Install Dependencies via System Package Manager

#### Ubuntu/Debian

```bash
sudo apt-get update
sudo apt-get install build-essential cmake git wget scons
sudo apt-get install libopenmpi-dev openmpi-bin
sudo apt-get install libhdf5-openmpi-dev
sudo apt-get install libboost-all-dev
sudo apt-get install libgsl-dev
sudo apt-get install libopenblas-dev  # Optimized BLAS for better performance
sudo apt-get install libblas-dev liblapack-dev
```

For PETSc and Exodus II, you may need to build from source (see Option 2).

#### CentOS/RHEL/Fedora

```bash
sudo yum install gcc gcc-c++ cmake git wget scons
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
brew install gcc cmake git wget scons
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

### Option 2: Build Dependencies from Source (Automated)

The repository includes a makefile to automatically download and build dependencies:

```bash
make all
```

This will:
- Download and build OpenMPI, HDF5, and Boost
- Install them in `$HOME/apollo/software`
- Configure build paths automatically

You can customize the installation by editing the makefile variables:
- `INSTALLDIR` - Where to install libraries (default: `$HOME/apollo/software`)
- `NUM_JOBS` - Parallel build jobs (default: 8)
- `CC`, `CXX` - Compiler choice

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

The executable will be created in `build-opt/` or `build-debug/` directory.

**Performance-Optimized Build:**

The optimized build (`build-opt`) includes Phase 1 performance optimizations that can provide **10-60x speedup**:

- **OpenMP threading**: Parallel element loop processing across CPU cores
- **Optimized BLAS**: Hardware-accelerated matrix operations (if OpenBLAS/CBLAS available)
- **Compiler optimizations**: `-O3 -march=native -ffast-math -ftree-vectorize`
- **Auto-vectorization**: SIMD instructions (AVX, AVX2, AVX-512)

To maximize performance:
1. Install OpenBLAS: `sudo apt-get install libopenblas-dev`
2. Set thread count: `export OMP_NUM_THREADS=8` (adjust to your CPU core count)
3. Build with: `scons build-opt`

**Note:** If dependencies are installed in non-standard locations, you may need to set environment variables:

```bash
export PATH=$HOME/apollo/software/bin:$PATH
export LD_LIBRARY_PATH=$HOME/apollo/software/lib:$LD_LIBRARY_PATH
export PETSC_DIR=$HOME/apollo/software/petsc
```

## Usage

### Basic Execution

Run simulations using MPI:

```bash
mpirun -np <number_of_processes> ./dg <input_file>
```

### Example Simulations

The `examples/unstructuredDG/` directory contains test cases for each physics module:

**Maxwell (Electromagnetics):**
```bash
mpirun -np 4 ./dg examples/unstructuredDG/maxwell/circularPulse/input.dat
```

**Euler (Gas Dynamics):**
```bash
mpirun -np 4 ./dg examples/unstructuredDG/euler/isentropicVortex/input.dat
```

**Multifluid (Plasma):**
```bash
mpirun -np 4 ./dg examples/unstructuredDG/multifluid/rmf_frc/input.dat
```

**Advection:**
```bash
mpirun -np 4 ./dg examples/unstructuredDG/advection/input.dat
```

## Example Problems

The repository includes various benchmark and application problems:

- **Maxwell:** Circular pulse, transverse magnetic modes, RMF antenna modeling for FRC
- **Euler:** Isentropic vortex, forward/backward facing steps, scramjet inlet, explosion tests
- **Multifluid:** Collisional-radiative models, RMF-driven FRC plasmas

Visualization scripts (Python) are provided in the `scripts/` directory for post-processing results.

## Performance

Apollo includes Phase 1 performance optimizations for high-performance computing:

### Implemented Optimizations

**1. OpenMP Thread Parallelism (4-8x speedup)**
- Element loop parallelized across CPU cores
- Shared-memory parallelism within each MPI rank
- Set threads: `export OMP_NUM_THREADS=<cores>` (e.g., 8 for an 8-core CPU)

**2. Optimized BLAS (2-5x speedup)**
- Custom matrix-vector operations replaced with hardware-accelerated BLAS
- Utilizes AVX/AVX2/AVX-512 SIMD instructions
- Requires OpenBLAS or similar library

**3. Compiler Optimizations (1.3-1.5x speedup)**
- Architecture-specific code generation (`-march=native`)
- Aggressive loop optimizations and vectorization
- Fast floating-point math

**Combined Expected Speedup: 10-60x** compared to unoptimized build

### Performance Tuning

**Thread Configuration:**
```bash
# For a 16-core CPU running 4 MPI ranks, use 4 threads per rank:
export OMP_NUM_THREADS=4
mpirun -np 4 ./build-opt/dg input.dat
```

**Hybrid MPI+OpenMP:**
- **Strong scaling**: Use more MPI ranks for small problems
- **Weak scaling**: Use OpenMP threads to utilize all cores per node
- **Recommended**: 1-4 MPI ranks per node, OpenMP for remaining cores

**Monitoring Performance:**
```bash
# Enable PETSc performance logging
export PETSC_OPTIONS="-log_view"
./build-opt/dg input.dat

# Profile with perf
perf stat -e cycles,instructions,cache-misses ./build-opt/dg input.dat
```

### Future Optimization Opportunities

Phase 2 and 3 optimizations (not yet implemented) could provide additional speedup:
- Communication/computation overlap (1.5-2x)
- GPU acceleration (10-50x for large problems)
- Mixed precision arithmetic (1.5-2x)
- Cache-friendly element reordering (1.2-1.5x)

## About PETSc

Apollo uses PETSc (Portable, Extensible Toolkit for Scientific Computation) for parallel numerical computations. PETSc provides:

- Highly portable parallel framework (distributed memory, shared memory, hybrid MPI-GPU)
- Comprehensive suite of parallel linear/nonlinear solvers
- Scalability to millions/billions of unknowns
- Industry-standard MPI-based communication
- Extensive documentation and active development community

## Documentation

Documentation is currently in development. For now, refer to:
- Example input files in `examples/unstructuredDG/`
- Visualization scripts in `scripts/`
- Source code headers in `src/`

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Author

Created by Eder Sousa