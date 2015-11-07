##
# Add PETSC libraries to build warpMConstructionEnvironment
##

Import('warpMConstructionEnv')


#Add specified library paths to search paths for libs and headers
if warpMConstructionEnv['petsc_base'] == '':
	warpMConstructionEnv['petsc_base']  = warpMConstructionEnv['$HOME/software'] + '/petsc-full'
	
if warpMConstructionEnv['petsc_base'] != '':
	#add base directory to path
	warpMConstructionEnv.AppendUnique(LIBPATH=(warpMConstructionEnv['petsc_base']+'/lib'))
	warpMConstructionEnv.AppendUnique(CPPPATH=(warpMConstructionEnv['petsc_base']+'/include'))


#configure libraries using some autoconf like functionality of scons
conf = Configure(warpMConstructionEnv)
print "PETSc_base set as ", warpMConstructionEnv['petsc_base']

## check petsc headers and libraries
if not conf.CheckCXXHeader('petscdmplex.h'):
        print "Header 'petscdmplex.h' was not found compile this program!"
	Exit(1)

if not conf.CheckCXXHeader('petsc.h'):
        print "Header 'petsc.h' was not found!"
        Exit(1)

#if not conf.CheckLib('MOAB'):
#        print "You need MOAB library to compile this program!"
#	Exit(1)

if not conf.CheckLib('petsc'):
        print "You need petsc library to compile this program!"
        Exit(1)

if not conf.CheckLib('exodus'):
        print "You need exodus library to compile this program!"
        Exit(1)

if not conf.CheckLib('MOAB'):
        print "You need MOAB library to compile this program!"
        Exit(1)        

#if not conf.CheckLib('exoIIv2for'):
#        print "You need exodus library to compile this program!"
#        Exit(1)

#if not conf.CheckLib('iMesh'):
#        print "You need exodus library to compile this program!"
#        Exit(1)

warpMConstructionEnv = conf.Finish() # replace the environment with the one modified by Configure's auto-conf actions
