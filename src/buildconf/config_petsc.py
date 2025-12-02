##
# Add PETSC libraries to build warpMConstructionEnvironment
##

Import('warpMConstructionEnv')
import os
import platform


#Add specified library paths to search paths for libs and headers
if warpMConstructionEnv['petsc_base'] == '':
	# Use environment variable if set, otherwise let SCons use system paths
	if 'PETSC_DIR' in os.environ:
		warpMConstructionEnv['petsc_base'] = os.environ['PETSC_DIR']
	# Don't set a default - let the system find PETSc

if warpMConstructionEnv['petsc_base'] != '':
	#add base directory to path
	warpMConstructionEnv.AppendUnique(LIBPATH=(warpMConstructionEnv['petsc_base']+'/lib'))
	warpMConstructionEnv.AppendUnique(CPPPATH=(warpMConstructionEnv['petsc_base']+'/include'))


#configure libraries using some autoconf like functionality of scons
conf = Configure(warpMConstructionEnv)
print("PETSc_base set as", warpMConstructionEnv['petsc_base'])

## check petsc headers and libraries
if not conf.CheckCXXHeader('petscdmplex.h'):
        print("Header 'petscdmplex.h' was not found compile this program!")
        Exit(1)

if not conf.CheckCXXHeader('petsc.h'):
        print("Header 'petsc.h' was not found!")
        Exit(1)

#if not conf.CheckLib('MOAB'):
#        print("You need MOAB library to compile this program!")
#        Exit(1)

if not conf.CheckLib('petsc'):
        print("You need petsc library to compile this program!")
        Exit(1)

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
