# -*- coding: utf-8 -*-
"""
Created on Mon Set 15 14:58:45 2015

@author: sousae
"""

#from scipy.interpolate import griddata
import wxunsdgdata2
from numpy import *
from optparse import OptionParser
from joblib import Parallel, delayed
from time import time
import multiprocessing
import pyvtk

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
                  help = 'Polynomial interpolation function order.',
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


def generateVTUfile(n):
	startT = time()
	spOrd = int(options.spatialOrder)
	filename = options.inputFile

	dh = wxunsdgdata2.WxVisData(filename,n)
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

	
	tris = range(3*dd.TotNumElements)
	
	tris = reshape(tris,(-1,3))
	
	grid = pyvtk.UnstructuredGrid(dd.gridPoints, triangle=tris)
	
	ptsData = pyvtk.PointData(pyvtk.Scalars(re, name="elec_rho"),
                           pyvtk.Scalars(vel_ex, name="elec_ux"),
                           pyvtk.Scalars(vel_ey, name="elec_vy"),
                           pyvtk.Scalars(vel_ez, name="elec_wz"),
                           pyvtk.Scalars(ee, name="elec_en"),
                           pyvtk.Scalars(ri, name="ion_rho"),
                           pyvtk.Scalars(vel_ix, name="ion_ux"),
                           pyvtk.Scalars(vel_iy, name="ion_vy"),
                           pyvtk.Scalars(vel_iz, name="ion_wz"),
                           pyvtk.Scalars(ei, name="ion_en"),
                           pyvtk.Scalars(Efieldx, name="Ex"),
                           pyvtk.Scalars(Efieldy, name="Ey"),
                           pyvtk.Scalars(Efieldz, name="Ez"),
                           pyvtk.Scalars(Bfieldx, name="Bx"),
                           pyvtk.Scalars(Bfieldy, name="By"),
                           pyvtk.Scalars(Bfieldz, name="Bz"))
	
	vtk = pyvtk.VtkData(grid, ptsData, 'Two-Fluid Plasma')
	
	#vtk = VtkData(\
    #    UnstructuredGrid(dd.gridPoints,
    #                     triangle=tris),
    #    PointData(Scalars(re,name='one')),
    #    PointData(Scalars(ee,name='two')),
    #    'Two-Fluid Plasma Unstructured Grid'
    #   )
	
	outfile = filename + '_2Fluid_' + str('%03d' % n)
	vtk.tofile(outfile,'binary')
	endT = time()
	totTime = endT - startT
	print("Frame "+str("%d" % n)+" COMPLETED in "+str("%3.2f" % totTime)+" secs.")

frame = int(options.frame)
stt = int(options.startFrame)
inputs = range(stt,frame+1)
num_cores = int(options.num_cores)
print("Generating plots using "+str("%d" % num_cores)+" processors.")
Parallel(n_jobs=num_cores)(delayed(generateVTUfile)(n) for n in inputs)
#generateVTUfile(0)
