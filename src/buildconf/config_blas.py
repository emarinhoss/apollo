##
# Add Blas libraries to build warpMConstructionEnvironment
##

Import('warpMConstructionEnv')

#Add specified library paths to search paths for libs and headers
if warpMConstructionEnv['blas_base'] == '':
        # Try common system locations
        warpMConstructionEnv['blas_base'] = '/usr/lib'

if warpMConstructionEnv['blas_base'] != '':
	#add base directory to path
        warpMConstructionEnv.AppendUnique(LIBPATH=(warpMConstructionEnv['blas_base']))
        # Also add x86_64-linux-gnu subdirectory common on Ubuntu/Debian
        warpMConstructionEnv.AppendUnique(LIBPATH=(warpMConstructionEnv['blas_base'] + '/x86_64-linux-gnu'))

#configure libraries using some autoconf like functionality of scons
conf = Configure(warpMConstructionEnv)
print "Blas_base set as ", warpMConstructionEnv['blas_base']

# Try to find BLAS library (try multiple common names)
blas_found = False

# Try OpenBLAS first (most common optimized BLAS)
if conf.CheckLib('openblas'):
    print "Found OpenBLAS library"
    blas_found = True
# Try standard BLAS + CBLAS
elif conf.CheckLib('cblas') or conf.CheckLib('blas'):
    print "Found standard BLAS library"
    if conf.CheckLib('cblas'):
        print "Found CBLAS library"
    blas_found = True
# Try combined BlasLapack (PETSc-provided)
elif conf.CheckLib('BlasLapack'):
    print "Found BlasLapack library"
    blas_found = True

if not blas_found:
    print "WARNING: No BLAS library found. Code will use slower fallback implementation."
    print "For best performance, install OpenBLAS: sudo apt-get install libopenblas-dev"
    # Don't exit - allow compilation without BLAS, will use fallback

# Check for gfortran (needed by some BLAS implementations)
if not conf.CheckLib('gfortran'):
    print "Note: gfortran not found, may be needed for some BLAS implementations"
    # Don't exit - not always needed

warpMConstructionEnv = conf.Finish() # replace the environment with the one modified by Configure's auto-conf actions
