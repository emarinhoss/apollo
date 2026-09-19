# -*- coding: utf-8 -*-
"""
Created on Mon Nov 24 15:42:19 2014

@author: sousae
"""
import os
import vtk
from numpy import array, zeros, append, reshape
from vtk.util.numpy_support import vtk_to_numpy

class WxDGArray:
    r"""WxDGArray(fname, frame, spOrd) -> WxDGArray

    WxArray objects behave like python arrays and hence support the
    python slice syntax. Note that only on array access is any data
    read from the HDF5 file."""

    def __init__(self, fname, spOrd):
        self.filename = fname
        self.order = spOrd
        # get pointer to data from open file handle
        #fq = arrayName
        #dps = 'self.wxdata.fh.root.%s' % comboSolver
        #dps = dps + '.__getattr__("%s")' % fq
        #self.dp = eval(dps)
        #self.dp.flavor = flavor
        #self.name = self.dp._v_name
        #self.fullShape = self.dp.shape
        #self.numComponents = self.fullShape[-1]
        #self.shape = self.fullShape[:-1]
        #self.onGrid = self.dp._v_attrs.vsMesh
        
        reader = vtk.vtkXMLUnstructuredGridReader()
        # tell the reader what is the filename to be read
        reader.SetFileName(self.filename)
        reader.Update()
        
        data = reader.GetOutput()
        self.nodesPerElem = (self.order+1)*(self.order+2)//2 # number of nodes per Element
        self.TotNumElements = 0 # Total number of Elements
        self.NumArrays  = data.GetCellData().GetNumberOfArrays() # Total number of arrays
        self.NumComp    = (self.NumArrays-1)//self.nodesPerElem # number of components per node
        #self.gridPoints = zeros((self.TotNumElements*self.nodesPerElem,3))
        #self.variables  = zeros((self.TotNumElements*self.nodesPerElem,self.NumComp))
        
        if(self.order==1):
            r = array([-1.,1.,-1.])
            s = array([-1.,-1.,1.])
            enum = 1
            connect = array([[0,1,2]])
        elif(self.order==2):
            r = array([-1.,0,1.,-1.,0,-1.])
            s = array([-1.,-1.,-1.,0,0.,1.])
            enum = 4
            connect = array([[0,1,3],[1,2,4],[1,4,3],[3,4,5]])
        elif(self.order==3):
            r = array([-1.,-0.447213595499958,0.447213595499958,1,-1,-0.333333333333333,0.447213595499958,-1.,-0.447213595499958,-1.])
            s = array([-1.,-1.,-1.,-1.,-0.447213595499958,-0.333333333333333,-0.447213595499958,0.447213595499958,0.447213595499958,1.])
            enum = 9
            connect = array([[0,5,4],[0,1,5],[1,2,5],[2,6,5],[2,3,6],[4,5,7],[5,8,7],[5,6,8],[7,8,9]])
        
        #r, s = elemtNodalPoints(self.order)
        coords = array([])
        varbls = array([])
        pts = zeros((self.nodesPerElem,2))
        
        for k in range(data.GetNumberOfCells()):
            cid=data.GetCell(k)
            
            p1 = data.GetPoint(cid.GetPointId(0))
            p2 = data.GetPoint(cid.GetPointId(1))
            p3 = data.GetPoint(cid.GetPointId(2))
            
            for np in range(len(r)):
                pts[np,0] = 0.5*(-p1[0]*(r[np]+s[np]) + p2[0]*(1.+r[np]) + p3[0]*(1.+ s[np]))
                pts[np,1] = 0.5*(-p1[1]*(r[np]+s[np]) + p2[1]*(1.+r[np]) + p3[1]*(1.+ s[np]))

            #tri = Delaunay(pts)
            #connect = tri.simplices.copy()
            self.TotNumElements += enum # Total number of Elements
            for kk in range(enum):
                for mm in  range(3):
                    coords = append(coords,pts[connect[kk,mm],0])
                    coords = append(coords,pts[connect[kk,mm],1])
                    coords = append(coords,0.0)
                    for cmps in range(self.NumComp):
                        value = vtk_to_numpy(data.GetCellData().GetArray(self.NumComp*connect[kk,mm]+cmps+1))
                        varbls = append(varbls,value[k])
            
        self.gridPoints = reshape(coords,(-1,3))
        self.variables = reshape(varbls,(-1,self.NumComp))
        
        self.res = self.variables

class WxVisData:
    r"""WxData(base : string, frm : int, flavor : string) -> WxData

    Provides an interface to read data from a WarpX hyperbolic solver
    simulation with base name ``base`` and  ``frm``.  Optionally
    the array ``flavor`` can be specified to select the kind of array
    to use (one of numpy or numeric).
    """
    
    def __init__(self, base, frm, flavor='numpy'):
        self.base = base
        self.frame = frm
        self.flavor = flavor
        fn = base + "_%d.vtu" % frm

        # ensure file exist
        if not os.path.exists(fn):
            raise Exception("WxData::__init__ : Dump %d of run %s not exist" % (frm, base))
            
        self.fname = fn
        # read in simulation time
        # self.time = float(self.fh.root.timeData._v_attrs.time)

    def close(self):
        r"""close() -> None

        Closes the file
        """
        self.fh.close()
        
    def readDG(self, spOrd):
        r"""readDG(name : string, spOrd : int) -> WxDGArray

        Read an array from the output file for a discontinuous
        galerkin simulation with ``meqn`` number of equations.  If
        ``comboSolver`` is specified it should be the name of the top
        comboSolver in the simulation.
        """

        return WxDGArray(self.fname, spOrd)
