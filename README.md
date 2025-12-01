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

- C++ compiler (with C++11 support)
- MPI library (OpenMPI, MPICH, or equivalent)
- PETSc library (version 3.6.0 or later recommended)
- CMake (for building)

## Installation

### Using the Makefile

The repository includes a makefile to install dependencies:

```bash
make all
```

This will download and build the required libraries (MPI, HDF5, Boost) in `$HOME/apollo/software`.

### Building Apollo

Once dependencies are installed:

```bash
mkdir build
cd build
cmake .. -DPETSC_WITH_MPI=ON
make
```

This will create an executable file called `dg` in the build directory.

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