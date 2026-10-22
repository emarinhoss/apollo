##
# Add Boost libraries to build warpMConstructionEnvironment
##

Import('warpMConstructionEnv')
import os
import platform
import sys


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

# Explicit directories, if given, win over anything derived from boost_base.
# SConstruct advertises both; before this they were accepted and ignored.
if warpMConstructionEnv['boost_includedir'] != '':
        warpMConstructionEnv.AppendUnique(CPPPATH=[warpMConstructionEnv['boost_includedir']])
if warpMConstructionEnv['boost_libdir'] != '':
        warpMConstructionEnv.AppendUnique(LIBPATH=[warpMConstructionEnv['boost_libdir']])

# Homebrew keeps Boost outside the default include path.
if sys.platform == 'darwin' and warpMConstructionEnv['boost_base'] == '':
        for prefix in ('/opt/homebrew', '/usr/local'):
                incDir = os.path.join(prefix, 'include')
                if os.path.isfile(os.path.join(incDir, 'boost', 'format.hpp')):
                        warpMConstructionEnv.AppendUnique(CPPPATH=[incDir])
                        warpMConstructionEnv.AppendUnique(LIBPATH=[os.path.join(prefix, 'lib')])
                        break


#configure libraries using some autoconf like functionality of scons
conf = Configure(warpMConstructionEnv)
print("BOOST_base set as", warpMConstructionEnv['boost_base'])

## check boost headers and libraries.
# Boost is header-only as far as Apollo is concerned, but it is not optional:
# hyperapps/multifluid/myIO.h includes <boost/format.hpp> and MathX.h and the
# grid code use uBLAS and Boost.Graph, so the build fails without it. Probe the
# header the multifluid module actually reaches for first, so the diagnostic
# names something the user can act on.
for header in ('boost/format.hpp',
               'boost/numeric/ublas/matrix.hpp',
               'boost/graph/adjacency_list.hpp',
               'boost/lexical_cast.hpp'):
        if conf.CheckCXXHeader(header):
                continue
        print("")
        print("ERROR: the Boost header <%s> was not found." % header)
        print("")
        print("  hyperapps/multifluid requires Boost (format, uBLAS and Graph).")
        print("    Debian/Ubuntu  sudo apt-get install libboost-dev")
        print("    RHEL/Fedora    sudo dnf install boost-devel")
        print("    macOS          brew install boost")
        print("")
        print("  Or point the build at an existing install:")
        print("    scons build-opt boost_base=/path/to/boost")
        print("")
        Exit(1)

warpMConstructionEnv = conf.Finish() # replace the environment with the one modified by Configure's auto-conf actions
