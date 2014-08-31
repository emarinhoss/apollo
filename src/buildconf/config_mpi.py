##
# Find MPI library
##



Import('warpMConstructionEnv')


#Add specified library paths to search paths for libs and headers

if warpMConstructionEnv['mpi_base'] != '':
	#add base directory to path so that scons can call correct mpicc
	warpMConstructionEnv.PrependENVPath('PATH', warpMConstructionEnv['mpi_base']+'/bin')
	warpMConstructionEnv.AppendUnique(LIBPATH=(warpMConstructionEnv['mpi_base']+'/lib'))
	warpMConstructionEnv.AppendUnique(CPPPATH=(warpMConstructionEnv['mpi_base']+'/include'))
#testing
import os
output = os.popen('which mpicc').read()
print "Which mpicc output", output
print "mpi_base set as ", warpMConstructionEnv['mpi_base']
##print "path is ", warpMConstructionEnv['ENV']['PATH']

#configure libraries using some autoconf like functionality of scons
conf = Configure(warpMConstructionEnv)

# the only configuration change required for MPI is to put in place the MPI compiler wrappers, then verify they work
conf.env['CXX'] = 'mpicxx'
conf.env['CC'] = 'mpicc'

try:
	assert(conf.CheckCHeader( 'mpi.h'))
	assert(conf.CheckCXXHeader( 'mpi.h'))
except:
	print 'MPI header mpi.h not found'
	raise

warpMConstructionEnv = conf.Finish() # replace the environment with the one modified by Configure's auto-conf actions
