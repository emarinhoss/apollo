#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 12 11:23:11 2018

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
                  help = 'Spatial order of the polynomial interpolation function.',
                  default = 1)
parser.add_option('-s', '--start', action = 'store',
                  dest = 'startFrame',
                  help = 'First frame to start plotting.',
                  default = 0)
parser.add_option('-g', '--gamma', action = 'store',
                  dest = 'gas_gamma',
                  help = 'Ratio of specific heats.',
                  default = 1.4)
parser.add_option('-n', '--numProcs', action = 'store',
                  dest = 'num_cores',
                  help = 'Number of cores to use.',
                  default = multiprocessing.cpu_count())

(options, args) = parser.parse_args()


frame = int(options.frame)
stt = int(options.startFrame)
gm = double(options.gas_gamma)

def generateVTUfile(n):
	spOrd = int(options.spatialOrder)
	filename = options.inputFile

	dh = wxunsdgdata2.WxVisData(filename,n)
	dd = dh.readDG(spOrd)

	rho = dd.variables[:,0]
	rhou = dd.variables[:,1]
	rhov = dd.variables[:,2]
	rhow = dd.variables[:,3]
	e = dd.variables[:,4]

	u = rhou/rho
	v = rhov/rho
	w = rhow/rho

	p = (gm-1.)*(e - 0.5*rho*(u*u+v*v+w*w))
	
	# element Type
	tris = range(3*dd.TotNumElements)
	
	tris = reshape(tris,(-1,3))
	
	grid = pyvtk.UnstructuredGrid(dd.gridPoints, triangle=tris)
	
	ptsData = pyvtk.PointData(pyvtk.Scalars(rho, name="density"),
                           pyvtk.Scalars(u, name="x_velovity"),
                           pyvtk.Scalars(v, name="y_velovity"),
                           pyvtk.Scalars(w, name="z_velovity"),
                           pyvtk.Scalars(p, name="pressure"))
	
	vtk = pyvtk.VtkData(grid, ptsData, 'Euler equations')

	outfile = filename + '_Euler_' + str('%03d' % n)
	vtk.tofile(outfile,'binary')
	print("Frame "+str("%d" % n)+" COMPLETE.")

inputs = range(stt,frame+1)
num_cores = int(options.num_cores)
print("Generating plots using "+str("%d" % num_cores)+" processors.")
Parallel(n_jobs=num_cores)(delayed(generateVTUfile)(n) for n in inputs)