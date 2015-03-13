##
# Add Blas libraries to build warpMConstructionEnvironment
##

Import('warpMConstructionEnv')


#Add specified library paths to search paths for libs and headers
if warpMConstructionEnv['blas_base'] == '':
        warpMConstructionEnv['blas_base']  = warpMConstructionEnv['/usr']
	
if warpMConstructionEnv['blas_base'] != '':
	#add base directory to path
        warpMConstructionEnv.AppendUnique(LIBPATH=(warpMConstructionEnv['blas_base']+'/lib'))


#configure libraries using some autoconf like functionality of scons
conf = Configure(warpMConstructionEnv)
print "Blas_base set as ", warpMConstructionEnv['blas_base']

if not conf.CheckLib('blas'):
        print "You need Blas library to compile this program!"
        Exit(1)

if not conf.CheckLib('lapack'):
        print "You need Lapack library to compile this program!"
        Exit(1)

warpMConstructionEnv = conf.Finish() # replace the environment with the one modified by Configure's auto-conf actions
