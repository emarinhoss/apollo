##
# Add Boost libraries to build warpMConstructionEnvironment
##

Import('warpMConstructionEnv')
import os
import platform


#Add specified library paths to search paths for libs and headers
# Use system-detected paths instead of hardcoded absolute paths
if warpMConstructionEnv['boost_base'] == '':
        # Let SCons use default system library search paths
        pass

if warpMConstructionEnv['boost_base'] != '':
        #add base directory to path
        warpMConstructionEnv.AppendUnique(LIBPATH=(warpMConstructionEnv['boost_base']+'/lib'))
        warpMConstructionEnv.AppendUnique(CPPPATH=(warpMConstructionEnv['boost_base']+'/include'))

        # Add architecture-specific subdirectory if it exists (Ubuntu/Debian convention)
        arch_subdir = platform.machine() + '-linux-gnu'
        arch_lib_path = os.path.join(warpMConstructionEnv['boost_base'], 'lib', arch_subdir)
        if os.path.isdir(arch_lib_path):
            warpMConstructionEnv.AppendUnique(LIBPATH=(arch_lib_path))


#configure libraries using some autoconf like functionality of scons
conf = Configure(warpMConstructionEnv)
print "BOOST_base set as ", warpMConstructionEnv['boost_base']

## check petsc headers and libraries
if not conf.CheckCXXHeader('boost/numeric/ublas/matrix.hpp'):
        print "Header 'matrix.hpp' was not found compile this program"
        Exit(1)

warpMConstructionEnv = conf.Finish() # replace the environment with the one modified by Configure's auto-conf actions
