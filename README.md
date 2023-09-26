# Apollo

This repository contains a discontinuous Galerkin (DG) finite element method implementation for solving hyperbolic partial differential equations (PDEs) such as the Advection, Euler, and Maxwell equations. The code is designed to be flexible and extensible, and can be used to solve a wide range of problems in plasma physics, including nuclear fusion research and space propulsion.

**Parallelization**

The code is parallelized using PETSc (Portable, Extensible Toolkit for Scientific Computation). PETSc is a suite of data structures and routines for the scalable (parallel) solution of scientific applications modeled by partial differential equations. It employs the Message Passing Interface (MPI) standard for all message-passing communication.

**Usage**

To run the code in parallel, simply type the following command:

```
mpirun -np <number_of_processes> ./dg <input_file>
```

This will run the code on <number_of_processes> processors.

**Example**

To run the example problem `plasma_transport.dat` in parallel on 4 processors, type the following command:

```
mpirun -np 4 ./dg examples/plasma_transport/input.dat
```

This will create an output file called `output.dat` in the `examples/plasma_transport/` directory.

**Requirements**

* C++ compiler
* MPI library
* PETSc library

**Installation**

To install the code, simply clone the repository and run the following commands:

```
mkdir build
cd build
cmake .. -DPETSC_WITH_MPI=ON
make
```

This will create an executable file called `dg`.

**Documentation**

In-progress...

## What is PETSc?

PETSc (Portable, Extensible Toolkit for Scientific Computation) is a suite of data structures and routines for the scalable (parallel) solution of scientific applications modeled by partial differential equations. It employs the Message Passing Interface (MPI) standard for all message-passing communication. PETSc is the world's most widely used parallel numerical software library for partial differential equations and sparse matrix computations.

## Why use PETSc?

PETSc has a number of advantages over other parallelization libraries, including:

* It is highly portable and can be used on a wide range of platforms, including distributed memory machines, shared memory machines, and hybrid MPI-GPU systems.
* It provides a wide range of parallel algorithms and data structures for solving a variety of scientific problems.
* It is highly scalable and can be used to solve large-scale problems with millions or even billions of unknowns.
* It is well-documented and supported by a large community of users and developers.