
# This is where I will install the support libraries for warpxm
INSTALLDIR = $(HOME)/apollo/software

# This is where I will temporarily store the src files for the libraries while they are built
SRCDIR = $(INSTALLDIR)/src

# This is where I will store the tar files before extraction
STORAGEDIR = $(INSTALLDIR)

# The number of jobs to build our libraries with
NUM_JOBS = 8

# Set the base compilers
CC = gcc
CXX = g++

# Set the downloader
DWNLD = wget

# These are the library names, if updating to a more recent library you (in theory) just need to modify these names... and the download link
PETSCNAME=petsc-3.6.0

# Define the install directories
PETSCINSTALLDIR=$(INSTALLDIR)/petsc

# Define the source zip files
MPIZIP=$(STORAGEDIR)/$(MPINAME).tar.gz

# Define the source directories
MPISRCDIR=$(SRCDIR)/$(MPINAME)


# Define the MPI compilers that will be used to build hdf5, boost and metis
MPICC = /usr/bin/mpicc
MPICXX = /usr/bin/mpicxx

all: initialize builddep_mpi builddep_hdf5 builddep_boost builddep_uc

initialize: cleandeps
	mkdir -p $(SRCDIR)
	mkdir -p $(STORAGEDIR)

builddep_mpi:
	# Check if mpi is installed
ifeq ("$(wildcard $(MPIINSTALLDIR)/bin/mpicxx)","")
	
	# Get OpenMPI
ifeq ("$(wildcard $(MPIZIP))","")
	cd $(STORAGEDIR);$(DWNLD) http://www.open-mpi.org/software/ompi/v1.8/downloads/openmpi-1.8.3.tar.gz
endif
	
	# Extract mpi
	tar xfz $(MPIZIP) -C $(SRCDIR)
	
	# Install mpi 
	cd $(MPISRCDIR);CC=$(CC) CXX=$(CXX) ./configure --prefix=$(MPIINSTALLDIR) --enable-shared --disable-static;make all -j $(NUM_JOBS);make install
endif

cleandeps:
	rm -rf $(SRCDIR)
