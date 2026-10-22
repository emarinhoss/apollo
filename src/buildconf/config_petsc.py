##
# Add PETSC libraries to build warpMConstructionEnvironment
##

Import('warpMConstructionEnv')
import os
import platform


# Locate PETSc. Order of precedence:
#   1. petsc_includedir / petsc_libdir, if given explicitly
#   2. petsc_base (command line, then this_host_config.py)
#   3. $PETSC_DIR, honouring $PETSC_ARCH
#   4. the compiler's default search paths (a distro package such as
#      libpetsc-real-dev, or a PETSc built with --prefix)
if warpMConstructionEnv['petsc_base'] == '' and 'PETSC_DIR' in os.environ:
        warpMConstructionEnv['petsc_base'] = os.environ['PETSC_DIR']

petscBase = warpMConstructionEnv['petsc_base']
petscIncludeDirs = []
petscLibDirs = []

if petscBase != '':
        petscIncludeDirs.append(os.path.join(petscBase, 'include'))
        petscLibDirs.append(os.path.join(petscBase, 'lib'))

        # A PETSc configured without --prefix keeps the generated petscconf.h and
        # the built libraries under $PETSC_DIR/$PETSC_ARCH, not directly under
        # $PETSC_DIR. Without this, <petsc.h> is found but petscconf.h is not and
        # the build fails with a confusing missing-header error deep in PETSc.
        petscArch = os.environ.get('PETSC_ARCH', '')
        if petscArch and os.path.isdir(os.path.join(petscBase, petscArch)):
                petscIncludeDirs.append(os.path.join(petscBase, petscArch, 'include'))
                petscLibDirs.append(os.path.join(petscBase, petscArch, 'lib'))
                print("PETSC_ARCH set as", petscArch)

# Explicit include/lib directories win over anything derived from petsc_base.
# These two variables are advertised by SConstruct; before this they were
# accepted and then ignored.
if warpMConstructionEnv['petsc_includedir'] != '':
        petscIncludeDirs = [warpMConstructionEnv['petsc_includedir']]
if warpMConstructionEnv['petsc_libdir'] != '':
        petscLibDirs = [warpMConstructionEnv['petsc_libdir']]

for d in petscIncludeDirs:
        warpMConstructionEnv.AppendUnique(CPPPATH=[d])
for d in petscLibDirs:
        warpMConstructionEnv.AppendUnique(LIBPATH=[d])


#configure libraries using some autoconf like functionality of scons
conf = Configure(warpMConstructionEnv)
print("PETSc_base set as", warpMConstructionEnv['petsc_base'])

def petscNotFound(what):
        print("")
        print("ERROR: %s was not found." % what)
        print("")
        print("  Apollo requires PETSc built with MPI support.")
        print("    Debian/Ubuntu  sudo apt-get install libpetsc-real-dev")
        print("    macOS          brew install petsc")
        print("    from source    see the 'Install PETSc' section of README.md")
        print("")
        print("  Then point the build at it, in order of precedence:")
        print("    scons build-opt petsc_base=/path/to/petsc")
        print("    export PETSC_DIR=/path/to/petsc  [and PETSC_ARCH=... if not a")
        print("      --prefix install]")
        print("")
        print("  Searched include dirs: %s" % (petscIncludeDirs or ['<system defaults>']))
        print("  Searched library dirs: %s" % (petscLibDirs or ['<system defaults>']))
        print("")
        Exit(1)

## check petsc headers and libraries
if not conf.CheckCXXHeader('petsc.h'):
        petscNotFound("the PETSc header 'petsc.h'")

if not conf.CheckCXXHeader('petscdmplex.h'):
        petscNotFound("the PETSc header 'petscdmplex.h' (Apollo uses DMPlex for "
                      "unstructured meshes)")

#if not conf.CheckLib('MOAB'):
#        print("You need MOAB library to compile this program!")
#        Exit(1)

if not conf.CheckLib('petsc'):
        petscNotFound("the PETSc library (libpetsc)")

# Try to find exodus library (may be bundled with PETSc or separate)
exodus_found = False
if conf.CheckLib('exodus'):
    print("Found exodus library")
    exodus_found = True
elif conf.CheckLib('exoIIv2c'):
    print("Found exoIIv2c library (exodus variant)")
    exodus_found = True
elif conf.CheckLib('exodusII'):
    print("Found exodusII library")
    exodus_found = True

if not exodus_found:
    print("WARNING: exodus library not found. Mesh I/O may be limited.")
    print("To install exodus on macOS: brew install seacas")
    print("Or rebuild PETSc with: --download-exodusii=1")
    # Don't exit - allow compilation without exodus if mesh I/O not needed

# GSL is required
if not conf.CheckLib('gsl'):
        print("You need gsl library to compile this program!")
        Exit(1)

if not conf.CheckLib('gslcblas'):
        print("You need gslcblas library to compile this program!")
        Exit(1)

#if not conf.CheckLib('MOAB'):
#        print("You need MOAB library to compile this program!")
#        Exit(1)

#if not conf.CheckLib('exoIIv2for'):
#        print("You need exodus library to compile this program!")
#        Exit(1)

#if not conf.CheckLib('iMesh'):
#        print("You need exodus library to compile this program!")
#        Exit(1)

warpMConstructionEnv = conf.Finish() # replace the environment with the one modified by Configure's auto-conf actions
