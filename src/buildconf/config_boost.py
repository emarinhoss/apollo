##
# Add PETSC libraries to build warpMConstructionEnvironment
##

Import('warpMConstructionEnv')


#Add specified library paths to search paths for libs and headers
if warpMConstructionEnv['boost_base'] == '':
        warpMConstructionEnv['boost_base']  = warpMConstructionEnv['/usr']

if warpMConstructionEnv['boost_base'] != '':
        #add base directory to path
        warpMConstructionEnv.AppendUnique(LIBPATH=(warpMConstructionEnv['boost_base']+'/lib/x86_64-linux-gnu'))
        warpMConstructionEnv.AppendUnique(CPPPATH=(warpMConstructionEnv['boost_base']+'/include'))


#configure libraries using some autoconf like functionality of scons
conf = Configure(warpMConstructionEnv)
print "BOOST_base set as ", warpMConstructionEnv['boost_base']

## check petsc headers and libraries
if not conf.CheckCXXHeader('boost/numeric/ublas/matrix.hpp'):
        print "Header 'matrix.hpp' was not found compile this program"
        Exit(1)

warpMConstructionEnv = conf.Finish() # replace the environment with the one modified by Configure's auto-conf actions
