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

	Ex = dd.variables[:,0]
	Ey = dd.variables[:,1]
	Ez = dd.variables[:,2]

	Bx = dd.variables[:,3]
	By = dd.variables[:,4]
	Bz = dd.variables[:,5]
	
	phi = dd.variables[:,6]
	psi = dd.variables[:,6]
	
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
	
	ug.point_data.scalars = Ex
	ug.point_data.scalars.name = 'Ex'
	ug.point_data.add_array(Ey)
	ug.point_data.get_array(1).name = 'Ey'
	ug.point_data.add_array(Ey)
	ug.point_data.get_array(2).name = 'Ez'
	
	ug.point_data.add_array(Bx)
	ug.point_data.get_array(3).name = 'Bx'
	ug.point_data.add_array(By)
	ug.point_data.get_array(4).name = 'By'
	ug.point_data.add_array(Bz)
	ug.point_data.get_array(5).name = 'Bz'
	
	ug.point_data.add_array(phi)
	ug.point_data.get_array(6).name = 'phi'
	ug.point_data.add_array(psi)
	ug.point_data.get_array(7).name = 'psi'
	

	outfile = filename + '_Maxwell_' + str('%03d' % n)  + '.vtu'
	save_xml(ug, outfile)
	print "Frame "+str("%d" % n)+" COMPLETE."

inputs = range(stt,frame+1)
num_cores = int(options.num_cores)
print "Generating plots using "+str("%d" % num_cores)+" processors."
Parallel(n_jobs=num_cores)(delayed(generateVTUfile)(n) for n in inputs)
#generateVTUfile(0)
