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
from time import time
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
	startT = time()
	spOrd = int(options.spatialOrder)
	filename = options.inputFile

	dh = wxunsdgdata.WxVisData(filename,n)
	dd = dh.readDG(spOrd)

	re = dd.variables[:,0]
	vel_ex = dd.variables[:,1]/re
	vel_ey = dd.variables[:,2]/re
	vel_ez = dd.variables[:,3]/re
	ee = dd.variables[:,4]

	ri = dd.variables[:,5]
	vel_ix = dd.variables[:,6]/ri
	vel_iy = dd.variables[:,7]/ri
	vel_iz = dd.variables[:,8]/ri
	ei = dd.variables[:,9]

	Efieldx = dd.variables[:,10]
	Efieldy = dd.variables[:,11]
	Efieldz = dd.variables[:,12]

	Bfieldx = dd.variables[:,13]
	Bfieldy = dd.variables[:,14]
	Bfieldz = dd.variables[:,15]

	
	# number of nodes per element
	nodesP = (spOrd+1)*(spOrd+2)/2
	# element Type
	elem_type = tvtk.Triangle().cell_type
	tris = range(3*dd.TotNumElements)
	
	tris = reshape(tris,(-1,3))

	ug = tvtk.UnstructuredGrid(points=dd.gridPoints)
	ug.set_cells(elem_type, tris)

	ug.point_data.scalars = re
	ug.point_data.scalars.name = 'elec_rho'
	ug.point_data.add_array(ee)
	ug.point_data.get_array(1).name = 'elec_en'

	ug.point_data.add_array(ri)
	ug.point_data.get_array(2).name = 'ion_rho'
	ug.point_data.add_array(ei)
	ug.point_data.get_array(3).name = 'ion_en'

	ug.point_data.add_array(vel_ex)
	ug.point_data.get_array(4).name = 'elec_ux'
	ug.point_data.add_array(vel_ey)
	ug.point_data.get_array(5).name = 'elec_vy'
	ug.point_data.add_array(vel_ez)
	ug.point_data.get_array(6).name = 'elec_wz'

	ug.point_data.add_array(vel_ix)
	ug.point_data.get_array(7).name = 'ion_ux'
	ug.point_data.add_array(vel_iy)
	ug.point_data.get_array(8).name = 'ion_vy'
	ug.point_data.add_array(vel_iz)
	ug.point_data.get_array(9).name = 'ion_wz'

	ug.point_data.add_array(Efieldx)
	ug.point_data.get_array(10).name = 'Ex'
	ug.point_data.add_array(Efieldy)
	ug.point_data.get_array(11).name = 'Ey'
	ug.point_data.add_array(Efieldz)
	ug.point_data.get_array(12).name = 'Ez'

	ug.point_data.add_array(Bfieldx)
	ug.point_data.get_array(13).name = 'Bx'
	ug.point_data.add_array(Bfieldy)
	ug.point_data.get_array(14).name = 'By'
	ug.point_data.add_array(Bfieldz)
	ug.point_data.get_array(15).name = 'Bz'

	outfile = filename + '_2Fluid_' + str('%03d' % n)  + '.vtu'
	save_xml(ug, outfile)
	endT = time()
	totTime = endT - startT
	print "Frame "+str("%d" % n)+" COMPLETED in "+str("%3.2f" % totTime)+" secs."

inputs = range(stt,frame+1)
num_cores = int(options.num_cores)
print "Generating plots using "+str("%d" % num_cores)+" processors."
Parallel(n_jobs=num_cores)(delayed(generateVTUfile)(n) for n in inputs)
#generateVTUfile(0)
