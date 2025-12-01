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
sudo apt-get install build-essential cmake git wget
sudo apt-get install libopenmpi-dev openmpi-bin
sudo apt-get install libhdf5-openmpi-dev
sudo apt-get install libboost-all-dev
sudo apt-get install libgsl-dev
sudo apt-get install libblas-dev liblapack-dev
```

For PETSc and Exodus II, you may need to build from source (see Option 2).

#### CentOS/RHEL/Fedora

```bash
sudo yum install gcc gcc-c++ cmake git wget
sudo yum install openmpi openmpi-devel
sudo yum install hdf5-openmpi hdf5-openmpi-devel
sudo yum install boost boost-devel
sudo yum install gsl gsl-devel
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
brew install gcc cmake git wget
brew install open-mpi
brew install hdf5-mpi
brew install boost
brew install gsl
brew install openblas lapack
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

Once all dependencies are installed:

```bash
cd /path/to/apollo
mkdir build
cd build
cmake .. -DPETSC_WITH_MPI=ON
make -j8
```

This will create an executable file called `dg` in the build directory.

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