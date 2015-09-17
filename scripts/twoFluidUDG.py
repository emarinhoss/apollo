# -*- coding: utf-8 -*-
"""
Created on Mon Set 15 14:58:45 2015

@author: sousae
"""

#from scipy.interpolate import griddata
import wxunsdgdata
from numpy import *
from tvtk.api import tvtk
from optparse import OptionParser
from joblib import Parallel, delayed
import multiprocessing

# set command line options
parser = OptionParser()
parser.add_option('-i', '--input', action = 'store',
                  dest = 'inputFile',
                  help = 'Base name of simulation')
parser.add_option('-f', '--frame', action = 'store',
                  dest = 'frame',
                  help = 'Last frame number to plot',
                  default = 0)
parser.add_option('-o', '--Order', action = 'store',
                  dest = 'spatialOrder',
                  help = 'Spatial order of the polynomial interpolation function.',
                  default = 1)
parser.add_option('-s', '--start', action = 'store',
                  dest = 'startFrame',
                  help = 'First frame to start plotting.',
                  default = 0)
parser.add_option('-n', '--numProcs', action = 'store',
                  dest = 'num_cores',
                  help = 'Number of cores to use.',
                  default = multiprocessing.cpu_count())

(options, args) = parser.parse_args()


def save_xml(ug, file_name):
    """Shows how you can save the unstructured grid dataset to a VTK
    XML file."""
    w = tvtk.XMLUnstructuredGridWriter(input=ug, file_name=file_name)
    w.write()


frame = int(options.frame)
stt = int(options.startFrame)

def generateVTUfile(n):
	spOrd = int(options.spatialOrder)
	filename = options.inputFile

	dh = wxunsdgdata.WxVisData(filename,n)
	dd = dh.readDG(spOrd)

	vel_e = zeros((3*dd.TotNumElements,3))
	re = dd.variables[:,0]
	vel_e[:,0] = dd.variables[:,1]/re
	vel_e[:,1] = dd.variables[:,2]/re
	vel_e[:,2] = dd.variables[:,3]/re
	ee = dd.variables[:,4]

	vel_i = zeros((3*dd.TotNumElements,3))
	ri = dd.variables[:,5]
	vel_i[:,0] = dd.variables[:,6]/ri
	vel_i[:,1] = dd.variables[:,7]/ri
	vel_i[:,2] = dd.variables[:,8]/ri
	ei = dd.variables[:,9]

	Efield = zeros((3*dd.TotNumElements,3))
	Efield[:,0] = dd.variables[:,10]
	Efield[:,1] = dd.variables[:,11]
	Efield[:,2] = dd.variables[:,12]

	Bfield = zeros((3*dd.TotNumElements,3))
	Bfield[:,0] = dd.variables[:,13]
	Bfield[:,1] = dd.variables[:,14]
	Bfield[:,2] = dd.variables[:,15]

	
	# number of nodes per element
	nodesP = (spOrd+1)*(spOrd+2)/2
	# element Type
	elem_type = tvtk.Triangle().cell_type
	tris = zeros((dd.TotNumElements,3),'int')
	sk = 0
	
	for K in range(dd.TotNumElements):
		for pots in range(3):
			tris[K,pots] = sk
			sk += 1

	ug = tvtk.UnstructuredGrid(points=dd.gridPoints)
	ug.set_cells(elem_type, tris)
	ug.point_data.scalars = re
	ug.point_data.scalars.name = 'elec_rho'
	ug.point_data.add_array(ee)
	ug.point_data.get_array(1).name = 'elec_en'#,'ion_rho','ion_en']
	ug.point_data.add_array(ri)
	ug.point_data.get_array(2).name = 'ion_rho'
	ug.point_data.add_array(ei)
	ug.point_data.get_array(3).name = 'ion_en'
	ug.point_data.vectors = vel_e
	ug.point_data.vectors.name = 'elec_vel'
	ug.point_data.add_array(vel_i)
	ug.point_data.get_array(5).name = 'ion_vel'
	ug.point_data.add_array(Efield)
	ug.point_data.get_array(6).name = 'Efield'
	ug.point_data.add_array(Bfield)
	ug.point_data.get_array(7).name = 'Bfield'
	outfile = filename + '_2Fluid_' + str('%03d' % n)  + '.vtu'
	save_xml(ug, outfile)
	print "Frame "+str("%d" % n)+" COMPLETE."

inputs = range(stt,frame+1)
num_cores = int(options.num_cores)
print "Generating plots using "+str("%d" % num_cores)+" processors."
Parallel(n_jobs=num_cores)(delayed(generateVTUfile)(n) for n in inputs)
#generateVTUfile(0)