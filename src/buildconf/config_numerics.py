##
# Add numerical libraries to build warpMConstructionEnvironment
##

Import('warpMConstructionEnv')


#Add specified library paths to search paths for libs and headers

if warpMConstructionEnv['boost_numeric_bindings_base'] == '':
	warpMConstructionEnv['boost_numeric_bindings_base']  = warpMConstructionEnv['RESEARCH_SOFTWARE_BASE'] + '/boost-numeric-bindings'

warpMConstructionEnv.AppendUnique(CPPPATH=(warpMConstructionEnv['boost_numeric_bindings_base']))

#configure libraries using some autoconf like functionality of scons
conf = Configure(warpMConstructionEnv)

try:
	# assuming all of the binding libraries are present if the lapack.hpp bindings are present
	assert(conf.CheckCXXHeader('boost/numeric/bindings/lapack.hpp'))
except:
	print('BOOST library or header not found')
	raise


# now ensure linking to the BLAS and LAPACK libraries
if warpMConstructionEnv['useAppleOpenCL']:   #serves as an awkwardly named indicator that this is a OS X build
	conf.env.AppendUnique(FRAMEWORKS = 'Accelerate')
else:
	try:
		assert(conf.CheckLib('blas'))
	except:
		print('BLAS library not found')
		raise
	try:
		assert(conf.CheckLib('lapack'))
	except:
		print('LAPACK library not found')
		raise


warpMConstructionEnv = conf.Finish() # replace the environment with the one modified by Configure's auto-conf actions
