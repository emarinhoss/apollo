# Convenience wrapper around the real build, which is SCons under src/.
#
# This file used to be an unfinished dependency bootstrapper: `make all` deleted
# $(SRCDIR) and then failed on three targets that were never written
# (builddep_hdf5, builddep_boost, builddep_uc), and its one real target expanded
# MPINAME, MPIZIP and MPIINSTALLDIR, none of which were ever defined. Since every
# dependency is now available as a distribution package - see README.md - the
# bootstrapper is not worth finishing, and these targets are useful instead.
#
# Everything here just calls what you would call by hand; nothing is hidden.

SCONS       ?= scons
PYTHON      ?= python3
SRCDIR      := src
JOBS        ?=
SCONSFLAGS  ?=

# Passed straight through to SCons, e.g. `make opt ARGS="arch=x86-64-v2 openmp=yes"`
ARGS ?=

.PHONY: help all opt debug both test test-python test-examples clean distclean deps-ubuntu

help:
	@echo "Apollo - make targets"
	@echo ""
	@echo "  make opt              build the optimized solver -> $(SRCDIR)/build-opt/apollo"
	@echo "  make debug            build the debug solver     -> $(SRCDIR)/build-debug/apollo"
	@echo "  make both             build both variants"
	@echo ""
	@echo "  make test             run every test (python + solver examples)"
	@echo "  make test-python      python tooling tests only (no compiler needed)"
	@echo "  make test-examples    run the solver against its bundled examples"
	@echo ""
	@echo "  make clean            remove build products"
	@echo "  make distclean        also remove SCons caches and preprocessed decks"
	@echo "  make deps-ubuntu      print the apt command that installs everything"
	@echo ""
	@echo "Options are passed to SCons through ARGS:"
	@echo "  make opt ARGS=\"arch=x86-64-v2\"     portable binary for older nodes"
	@echo "  make opt ARGS=\"petsc_base=/opt/petsc\""
	@echo "  $(SCONS) -h                            (run in $(SRCDIR)/) full variable list"

all: opt

opt:
	cd $(SRCDIR) && $(SCONS) $(SCONSFLAGS) $(JOBS) build-opt $(ARGS)

debug:
	cd $(SRCDIR) && $(SCONS) $(SCONSFLAGS) $(JOBS) build-debug $(ARGS)

both: opt debug

test: test-python test-examples

test-python:
	$(PYTHON) -m unittest discover -s test -v

test-examples:
	test/run_examples.sh

clean:
	cd $(SRCDIR) && $(SCONS) -c build-opt build-debug || true
	rm -rf $(SRCDIR)/build-opt $(SRCDIR)/build-debug

distclean: clean
	rm -rf $(SRCDIR)/.sconf_temp $(SRCDIR)/.sconsign.dblite $(SRCDIR)/config.log
	find . -name '__pycache__' -type d -prune -exec rm -rf {} +
	find examples -name '*.inp' -o -name '*_temp1' -o -name '*_temp2' | xargs -r rm -f

deps-ubuntu:
	@echo "sudo apt-get install -y build-essential git python3-pip \\"
	@echo "    libopenmpi-dev openmpi-bin libpetsc-real-dev libhdf5-openmpi-dev \\"
	@echo "    libboost-dev libeigen3-dev libgsl-dev libopenblas-dev"
	@echo "pip3 install scons"
