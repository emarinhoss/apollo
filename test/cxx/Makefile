# Unit tests for pieces of Apollo that can be exercised without a mesh, a
# solver, or MPI ranks.
#
#     make -C test/cxx          build and run everything
#     make -C test/cxx clean
#
# These link against the objects SCons already built, so build the solver first
# (cd src && scons build-opt). PETSC_DIR is honoured the same way the SCons
# build honours it.

REPO      := $(abspath $(dir $(lastword $(MAKEFILE_LIST)))/../..)
SRC       := $(REPO)/src
BUILD     ?= $(SRC)/build-opt

# Plain assignment, not ?=: make predefines CXX as g++, so ?= would never fire.
# A command-line "make CXX=..." still wins over this.
CXX        = mpicxx
PETSC_DIR ?= /usr/lib/petsc

# Mirror the include set from src/SConscript, plus the dependency headers.
INCLUDES := -I$(SRC)/lib -I$(SRC)/lapack_lite -I$(SRC)/solvers -I$(SRC)/subsolvers \
            -I$(SRC)/hyper -I$(SRC)/hyperapps -I$(SRC)/xapollo \
            -I$(PETSC_DIR)/include

# Where distributions put the headers that are not on the default path.
INCLUDES += $(addprefix -I,$(wildcard /usr/include/hdf5/openmpi /usr/include/eigen3 \
                                      /usr/local/include/eigen3))

# -MMD -MP emits a .d file listing every header the test pulled in. Without it
# make would rebuild only when the .cc changed - useless here, since most of the
# code under test lives in headers.
CXXFLAGS   = -std=c++14 -g -O1 -Wall -MMD -MP
LDFLAGS    = -L$(PETSC_DIR)/lib
LDLIBS     = -lpetsc -lgsl -lgslcblas

# Input-deck parsing, logging and the object base: needed by every test that
# configures a real class from a real deck fragment.
CORE_OBJS := $(BUILD)/lib/wxcryptset.o $(BUILD)/lib/wxcryptsetlexer.o \
             $(BUILD)/lib/wxcrypt.o $(BUILD)/lib/wxobject.o \
             $(BUILD)/lib/wxexcept.o $(BUILD)/lib/wxlogger.o \
             $(BUILD)/lib/wxlogstream.o $(BUILD)/lib/wxlogstreamstrm.o

# Per-test extras. Kept minimal on purpose: a link error should point at the
# code under test rather than drag in the whole solver.
test_friction_sources_OBJS := $(BUILD)/hyper/wxhyperbolicsrc.o
test_rmf_antenna_OBJS      := $(BUILD)/hyper/wxhyperbolicsrc.o

# The boundary-condition classes define their methods in .cc files with explicit
# template instantiation, so the object is needed even though the class is
# reached through its header.
test_rmf_boundary_OBJS := $(BUILD)/solvers/wxgridbc.o \
                          $(BUILD)/solvers/wxstepper.o \
                          $(BUILD)/subsolvers/apsubsolver.o \
                          $(BUILD)/hyperapps/multifluid/aptwofluidsimplifiedrmfbc.o

TESTS := test_friction_sources test_rmf_boundary test_rmf_antenna
BINS  := $(addprefix $(CURDIR)/,$(TESTS))

.PHONY: all run clean check-build

all: run

check-build:
	@test -d "$(BUILD)" || { \
	  echo "error: no build at $(BUILD)"; \
	  echo "       build the solver first:  cd src && scons build-opt"; \
	  echo "       or point this at another variant:  make -C test/cxx BUILD=$(SRC)/build-debug"; \
	  exit 2; }

$(CURDIR)/%: %.cc | check-build
	$(CXX) $(CXXFLAGS) $(INCLUDES) $< $(CORE_OBJS) $($*_OBJS) -o $@ $(LDFLAGS) $(LDLIBS)

-include $(addsuffix .d,$(TESTS))

run: $(BINS)
	@rc=0; for t in $(BINS); do \
	  echo "=== $$(basename $$t)"; \
	  $$t || rc=1; \
	done; exit $$rc

clean:
	rm -f $(BINS) $(addsuffix .d,$(TESTS))
