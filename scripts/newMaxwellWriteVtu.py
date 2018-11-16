#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 12 11:37:55 2018

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
parser.add_option('-n', '--numProcs', action = 'store',
                  dest = 'num_cores',
                  help = 'Number of cores to use.',
                  default = multiprocessing.cpu_count())

(options, args) = parser.parse_args()

frame = int(options.frame)
stt = int(options.startFrame)

def generateVTUfile(n):
	spOrd = int(options.spatialOrder)
	filename = options.inputFile

	dh = wxunsdgdata2.WxVisData(filename,n)
	dd = dh.readDG(spOrd)

	Ex = dd.variables[:,0]
	Ey = dd.variables[:,1]
	Ez = dd.variables[:,2]

	Bx = dd.variables[:,3]
	By = dd.variables[:,4]
	Bz = dd.variables[:,5]
	
	phi = dd.variables[:,6]
	psi = dd.variables[:,7]
    
    # element Type
	tris = range(3*dd.TotNumElements)
	
	tris = reshape(tris,(-1,3))
	
	grid = pyvtk.UnstructuredGrid(dd.gridPoints, triangle=tris)
	
	ptsData = pyvtk.PointData(pyvtk.Scalars(Ex, name="Ex"),
                           pyvtk.Scalars(Ey, name="Ey"),
                           pyvtk.Scalars(Ez, name="Ez"),
                           pyvtk.Scalars(Bx, name="Bx"),
                           pyvtk.Scalars(By, name="By"),
                           pyvtk.Scalars(Bz, name="Bz"),
                           pyvtk.Scalars(phi, name="phi"),
                           pyvtk.Scalars(psi, name="psi"))
	
	vtk = pyvtk.VtkData(grid, ptsData, 'Maxwell equations')
	

	outfile = filename + '_Maxwell_' + str('%03d' % n)
	vtk.tofile(outfile,'binary')
	print("Frame "+str("%d" % n)+" COMPLETE.")

inputs = range(stt,frame+1)
num_cores = int(options.num_cores)
print("Generating plots using "+str("%d" % num_cores)+" processors.")
Parallel(n_jobs=num_cores)(delayed(generateVTUfile)(n) for n in inputs)