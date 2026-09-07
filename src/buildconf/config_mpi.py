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
output = os.popen('which mpicxx').read()
print("Which mpicxx output", output)
print("mpi_base set as", warpMConstructionEnv['mpi_base'])
##print("path is", warpMConstructionEnv['ENV']['PATH'])

#configure libraries using some autoconf like functionality of scons
conf = Configure(warpMConstructionEnv)

# the only configuration change required for MPI is to put in place the MPI compiler wrappers, then verify they work
conf.env['CXX'] = 'mpicxx'
conf.env['CC'] = 'mpicc'

if not (conf.CheckCHeader('mpi.h') and conf.CheckCXXHeader('mpi.h')):
	print("")
	print("ERROR: could not compile against <mpi.h> using the mpicc/mpicxx wrappers.")
	print("")
	print("  Apollo needs an MPI development installation. Install one with:")
	print("    Debian/Ubuntu  sudo apt-get install libopenmpi-dev openmpi-bin")
	print("    RHEL/Fedora    sudo dnf install openmpi-devel   (then: module load mpi)")
	print("    macOS          brew install open-mpi")
	print("")
	print("  If MPI is installed somewhere non-standard, point the build at it:")
	print("    scons build-opt mpi_base=/path/to/mpi")
	print("  mpi_base is currently '%s'." % warpMConstructionEnv['mpi_base'])
	print("  See src/this_host_config.py.")
	print("")
	Exit(1)

warpMConstructionEnv = conf.Finish() # replace the environment with the one modified by Configure's auto-conf actions
