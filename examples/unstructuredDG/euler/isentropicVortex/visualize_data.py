# -*- coding: utf-8 -*-
"""
Created on Tue Dec 16 14:02:45 2014

@author: sousae
"""

#from scipy.interpolate import griddata
import wxunsdgdata
from numpy import *
from tvtk.api import tvtk
from optparse import OptionParser

# set command line options
parser = OptionParser()
parser.add_option('-i', '--input', action = 'store',
                  dest = 'inputFile',
                  help = 'Base name of simulation')
parser.add_option('-f', '--frame', action = 'store',
                  dest = 'frame',
                  help = 'Last frame number to plot',
                  default = 0)
parser.add_option('-v', '--variable', action = 'store',
                  dest = 'variable',
                  help = 'Index of variable to plot, q(variable)',
                  default = 0)
parser.add_option('-s', '--spOrder', action = 'store',
                  dest = 'spatialOrder',
                  help = 'Spatial order of the polynomial interpolation function.',
                  default = 1)                  

(options, args) = parser.parse_args()



def save_xml(ug, file_name):
    """Shows how you can save the unstructured grid dataset to a VTK
    XML file."""
    w = tvtk.XMLUnstructuredGridWriter(input=ug, file_name=file_name)
    w.write()


# polynomial expansion order
spOrd = int(options.spatialOrder)
component = int(options.variable)
frame = int(options.frame)
filename = options.inputFile

for n in range(0,frame+1):
	dh = wxunsdgdata.WxVisData(filename,n)
	dd = dh.readDG(spOrd)
	var= dd.variables[:,component]
	
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
	ug.point_data.scalars = var
	ug.point_data.scalars.name = 'Q_var'
	outfile = filename + '_Comp_' + str(component) + '_' + str('%03d' % n)  + '.vtu'
	save_xml(ug, outfile)
