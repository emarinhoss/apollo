# -*- coding: utf-8 -*-
"""
Created on Tue Dec 16 14:02:45 2014

@author: sousae
"""

#from scipy.interpolate import griddata
import wxunsdgdata
from numpy import *
from tvtk.api import tvtk

def save_xml(ug, file_name):
    """Shows how you can save the unstructured grid dataset to a VTK
    XML file."""
    w = tvtk.XMLUnstructuredGridWriter(input=ug, file_name=file_name)
    w.write()

filename = 'advection'
# polynomial expansion order
spOrd = 1

dh = wxunsdgdata.WxVisData(filename,0)
dd = dh.readDG(spOrd)
var= dd.variables[:,0]

#xmin = min(dd.gridPoints[:,0])
#xmax = max(dd.gridPoints[:,0])
#ymin = min(dd.gridPoints[:,1])
#ymax = max(dd.gridPoints[:,1])
#ss = shape(dd.gridPoints)
#res = round(sqrt(ss[0]))
#res = 100
#points = 1.e-8*random.rand(ss[0], 2)
#grid = dd.gridPoints+points
#grid_x, grid_y = np.mgrid[xmin:xmax:res*1j, ymin:ymax:res*1j]
#grid_z2 = griddata(grid, var, (grid_x, grid_y), method='cubic')
#imshow(grid_z2.T, extent=(0,1,0,1), origin='lower')
#show()

# number of nodes per element
nodesP = (spOrd+1)*(spOrd+2)/2
# element Type
elem_type = tvtk.Triangle().cell_type
# 
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

save_xml(ug, 'file.vtu')