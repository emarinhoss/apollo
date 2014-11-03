# -*- coding: utf-8 -*-
"""
Created on Thu Oct 23 07:59:11 2014

This is a test program that reads data from the vtu files for DG data that is
stacked at the cell centers and distributes them on the correspondent nodes.


@author: sousae
"""

from vtk import *
from numpy import *
from vtk.util.numpy_support import vtk_to_numpy
from pyvisfile.vtk import ( 
    UnstructuredGrid, DataArray,
    AppendedDataXMLGenerator,
    VTK_VERTEX, VF_LIST_OF_VECTORS, VF_LIST_OF_COMPONENTS)


SPORD = 1
meqn  = 5
nodesPer = (SPORD+1)*(SPORD+2)/2
#
r = [[-1,1,-1],[-1,0,1,-1,0,-1],[-1,-0.447213595499958,0.447213595499958,1,-1,-0.333333333333333,0.447213595499958,-1,-0.447213595499958,-1]]
s = [[-1,-1,1],[-1,-1,-1,0,0,1],[-1,-1,-1,-1,-0.447213595499958,-0.333333333333333,-0.447213595499958,0.447213595499958,0.447213595499958,1]]

reader = vtk.vtkXMLUnstructuredGridReader()
reader.SetFileName("forwardFacingStep_0.vtu")
reader.Update()
#
data = reader.GetOutput()
#values = data.GetCellData().GetArray(1)
#VV = vtk_to_numpy(values)

x = zeros(data.GetNumberOfCells()*nodesPer)
y = zeros(data.GetNumberOfCells()*nodesPer)
rh= zeros(data.GetNumberOfCells()*nodesPer)
ru= zeros(data.GetNumberOfCells()*nodesPer)
rv= zeros(data.GetNumberOfCells()*nodesPer)
rw= zeros(data.GetNumberOfCells()*nodesPer)
en= zeros(data.GetNumberOfCells()*nodesPer)
sk = 0

for k in range(data.GetNumberOfCells()):
    cid=data.GetCell(k)
    
    p1 = data.GetPoint(cid.GetPointId(0))
    p2 = data.GetPoint(cid.GetPointId(1))
    p3 = data.GetPoint(cid.GetPointId(2))
    
    for np in range(len(r[SPORD-1])):
        x[sk] = 0.5*(-p1[0]*(r[SPORD-1][np]+s[SPORD-1][np]) + p2[0]*(1.+r[SPORD-1][np]) + p3[0]*(1.+ s[SPORD-1][np]))
        y[sk] = 0.5*(-p1[1]*(r[SPORD-1][np]+s[SPORD-1][np]) + p2[1]*(1.+r[SPORD-1][np]) + p3[1]*(1.+ s[SPORD-1][np]))
        rho = vtk_to_numpy(data.GetCellData().GetArray(meqn*np))
        rh[sk] = rho[k]
        rhou= vtk_to_numpy(data.GetCellData().GetArray(meqn*np+1))
        ru[sk] = rhou[k]
        rhov= vtk_to_numpy(data.GetCellData().GetArray(meqn*np+2))
        rv[sk] = rhov[k]
        rhow= vtk_to_numpy(data.GetCellData().GetArray(meqn*np+3))
        rw[sk] = rhow[k]
        ener= vtk_to_numpy(data.GetCellData().GetArray(meqn*np+4))
        en[sk] = ener[k]
        sk += 1
        

n = len(x)

points = zeros((n,3))
points[:,0] = x
points[:,1] = y


grid = UnstructuredGrid(
        (n, DataArray("points", points, vector_format=VF_LIST_OF_VECTORS)),
        cells=arange(n, dtype=uint32),
        cell_types=asarray([VTK_VERTEX] * n, dtype=uint8))

momen = zeros((n,3))
momen[:,0] = ru
momen[:,1] = rv
momen[:,2] = rw

data = [
        ("rho", rh),
        ("rhou", ru),
        ("rhov", rv),
        ("rhow", rw),
        ("energy", en),
]

file_name = "points.vts"
compressor = None

for name, field in data:
    grid.add_pointdata(DataArray(name, field,
        vector_format=VF_LIST_OF_COMPONENTS))
        
outf = open(file_name, "w")
AppendedDataXMLGenerator(compressor)(grid).write(outf)
outf.close()        