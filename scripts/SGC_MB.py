import os
import os.path
import sys
#import h5py
from numpy import *
import numpy as np
import random
import pdb
from time import clock, time
from copy import copy, deepcopy

sys.setrecursionlimit(1000000000) 

# This program takes an 'Abacus' File (.inp) from Cubit
# and converts it to an hdf5 file that only contains
# a 1, 2 or 3D matrix of nodepositions

# The following are classes that contain the Grid

# A Grid is made up of a group of domains

# Domains are fully contained structured Grids (with boundary ghost cells)

# This program functions by mapping out the connectivity of the Elements
# Then a random Element is 

class Node(object):

    # Nodes are the 'vertexes' of the Grid
    
    def __init__(self, ident, position):
        self.id = ident                 # Unique Identity of Node (global ID from Cubit)
        self.position = list(position)  # Physical position of Node
        self.indices = []               # Holds the Indices of the Node within the Grid
	self.neighbors = []             # Holds the neighboring cell originally found in the elements list
	self.geometry = []              # Holds the corresponding nodes associated with the cell 
	self.faces = []                 # Holds the corresponding faces associated with the cell, 
					# arranged in paired sets, i.e. each direction is presented 
					# as a set of two faces on opposite sides of the cell.
        self.domain = []
        self.celldata = []

    # Returns the Global ID of the Node
    def getID(self):
        return self.id

    # Returns the Local ID of the Node (-1 if Node is not in Domain)
    def getLocalID(self, domain_index = -1):
        if domain_index < 0:
            return -1
        else:
            if domain_index in self.owners:    
                return self.owners[domain_index][patch_index]
            else:
                return -1

    def getNeighbors(self):
        return self.neighbors

###########################################        
# ---------------------------------------------

    def findNumofSharedCells(self,externaldomain):

        neighbordata = self.neighbors
        neighbordomainlist = self.domain
        sharedcellslist = []

        for currentneighbor in range(len(neighbordata)):
            neighborcellid = neighbordata[currentneighbor]
            currentneighbordomain = neighbordomainlist[currentneighbor]
            domainmatches = (externaldomain==currentneighbordomain)
            if (domainmatches):
                sharedcellslist.append(neighborcellid)
#                print "domainmatches = " +str(domainmatches)
#                print "sharedcellslist @ nodeID ("+str(self.id)+ "= " +str(sharedcellslist)

        numofsharedcellsingivendomain= len(sharedcellslist)

        return numofsharedcellsingivendomain

# ---------------------------------------------

    def findSharedCellsForSingleDomain(self,externaldomain):

        neighbordata = self.neighbors
        neighbordomainlist = self.domain
        sharedcellslist = []

        for currentneighbor in range(len(neighbordata)):
            neighborcellid = neighbordata[currentneighbor]
            currentneighbordomain = neighbordomainlist[currentneighbor]
            domainmatches = (externaldomain==currentneighbordomain)
            if (domainmatches):
                sharedcellslist.append(neighborcellid)

        return sharedcellslist

# ---------------------------------------------

    def findSharedCellsGeometryForSingleDomain(self,externaldomain):

        neighbordata = self.neighbors
        neighbordomainlist = self.domain
        neighborgeometry = self.geometry
        sharedcellsgeometry = []

        for currentneighbor in range(len(neighbordata)):
#            neighborcellid = neighbordata[currentneighbor]
            neighborcellgeometry = neighborgeometry[currentneighbor]
            currentneighbordomain = neighbordomainlist[currentneighbor]
            domainmatches = (externaldomain==currentneighbordomain)
            if (domainmatches):
                sharedcellsgeometry.append(neighborcellgeometry)

        return sharedcellsgeometry

# ---------------------------------------------

    def findSharedCellsFaceListForSingleDomain(self,externaldomain):


        neighborfaces       = self.faces
        neighborcell_list   = self.neighbors
        neighbordomain_list = self.domain

        numofneighbors      = len(neighborcell_list)

        sharedcellsfacelist = []

        for currentneighbor in range(numofneighbors):
#            currentcellID         = neighborcell_list[currentneighbor]
            neighborcellfacelist  = neighborfaces[currentneighbor]
            currentneighbordomain = neighbordomain_list[currentneighbor]
            domainmatches = (externaldomain==currentneighbordomain)
            if (domainmatches):
                sharedcellsfacelist.append(neighborcellfacelist)

        return sharedcellsfacelist

# ---------------------------------------------

    def replaceFaceListForSpecifiedCell(self,externaldomainID,externalface_list,externalcellID):

        neighborcell_list = self.neighbors
        neighbordomain_list = self.domain
        numofneighbors = len(neighborcell_list)

        for externalface in externalface_list:
            checksum = 0
            for node in range(len(externalface)):
                checksum = checksum+externalface[node]
            if (checksum == 0):
                break

        if (checksum != 0):
            for currentneighbor in range(numofneighbors):
                currentcellID = neighborcell_list[currentneighbor]
                currentdomainID = neighbordomain_list[currentneighbor]
                domainmatches = (externaldomainID==currentdomainID)
                cellmatches = (externalcellID==currentcellID)
                if (domainmatches and cellmatches):
                    self.faces[currentneighbor]=externalface_list

# ---------------------------------------------
# ---------------------------------------------

    def replaceGeometryListForSpecifiedCell(self,externaldomainID,externalgeometry_list,externalcellID):

        neighborcell_list = self.neighbors
        neighbordomain_list = self.domain
        numofneighbors = len(neighborcell_list)

        for currentneighbor in range(numofneighbors):
            currentcellID = neighborcell_list[currentneighbor]
            currentdomainID = neighbordomain_list[currentneighbor]
            domainmatches = (externaldomainID==currentdomainID)
            cellmatches = (externalcellID==currentcellID)
            if (domainmatches and cellmatches):
                self.geometry[currentneighbor]=externalgeometry_list

# ---------------------------------------------

    def printSpecificCellForASpecificNode(self,nodes,externalcellID):

        neighborcell_list = self.neighbors
        numofneighbors = len(neighborcell_list)

#        cellmatches = False

        for currentneighbor in range(numofneighbors):
            currentcellID = neighborcell_list[currentneighbor]
            cellfound = (externalcellID==currentcellID)
            if (cellfound):
                break
        if (cellfound):
            facelist = self.faces[currentneighbor]
            xpos = facelist[0]
            xneg = facelist[1]

            n000 = nodes[xneg[0]]
            nx00 = nodes[xpos[0]]
            n0y0 = nodes[xneg[3]]
            n00z = nodes[xneg[1]]
            nxy0 = nodes[xpos[3]]
            n0yz = nodes[xneg[2]]
            nx0z = nodes[xpos[1]]
            nxyz = nodes[xpos[2]]

            n000pos = n000.position
            nx00pos = nx00.position
            n0y0pos = n0y0.position
            n00zpos = n00z.position
            nxy0pos = nxy0.position
            n0yzpos = n0yz.position
            nx0zpos = nx0z.position
            nxyzpos = nxyz.position



            print 'Given nnID = ' + str(self.id) + ', nodes for cell '+str(externalcellID)+':'
            print 'n000 is '+str(n000.id)+ ' = ('+str(round(n000pos[0],1))+', \t'+str(round(n000pos[1],1))+', \t'+str(round(n000pos[2],1))+')'
            print 'nx00 is ' +str(nx00.id)+ ' = ('+str(round(nx00pos[0],1))+', \t'+str(round(nx00pos[1],1))+', \t'+str(round(nx00pos[2],1))+')'
            print 'n0y0 is ' +str(n0y0.id)+ ' = ('+str(round(n0y0pos[0],1))+', \t'+str(round(n0y0pos[1],1))+', \t'+str(round(n0y0pos[2],1))+')'
            print 'n00z is ' +str(n00z.id)+ ' = ('+str(round(n00zpos[0],1))+', \t'+str(round(n00zpos[1],1))+', \t'+str(round(n00zpos[2],1))+')'
            print 'nxy0 is ' +str(nxy0.id)+ ' = ('+str(round(nxy0pos[0],1))+', \t'+str(round(nxy0pos[1],1))+', \t'+str(round(nxy0pos[2],1))+')'
            print 'n0yz is ' +str(n0yz.id)+ ' = ('+str(round(n0yzpos[0],1))+', \t'+str(round(n0yzpos[1],1))+', \t'+str(round(n0yzpos[2],1))+')'
            print 'nx0z is ' +str(nx0z.id)+ ' = ('+str(round(nx0zpos[0],1))+', \t'+str(round(nx0zpos[1],1))+', \t'+str(round(nx0zpos[2],1))+')'
            print 'nxyz is ' +str(nxyz.id)+ ' = ('+str(round(nxyzpos[0],1))+', \t'+str(round(nxyzpos[1],1))+', \t'+str(round(nxyzpos[2],1))+')'
        else:
            print 'for nodeID: '+str(self.id)+':  cellID = '+str(externalcellID)+' not found'
# ---------------------------------------------

    def calculateFaces(self, style, geometrylist):

        if (style == 'Hexahedron'):
            face1pindices,face1nindices,face2pindices,face2nindices,face3pindices,face3nindices = cellHexahedronFaceIndices()
            face1p = [geometrylist[i] for i in face1pindices]
            face1n = [geometrylist[i] for i in face1nindices]
            face2p = [geometrylist[i] for i in face2pindices]
            face2n = [geometrylist[i] for i in face2nindices]
            face3p = [geometrylist[i] for i in face3pindices]
            face3n = [geometrylist[i] for i in face3nindices]
            faces = (face1p,face1n,face2p,face2n,face3p,face3n)
        elif (style == 'Rectangle'):
            face1nindices = (0,1)
            face2nindices = (0,3)
            face1pindices = (3,2)
            face2pindices = (1,2)
            face1p = [geometrylist[i] for i in face1pindices]
            face1n = [geometrylist[i] for i in face1nindices]
            face2p = [geometrylist[i] for i in face2pindices]
            face2n = [geometrylist[i] for i in face2nindices]
            faces = (face1p,face1n,face2p,face2n)

        return faces
# ---------------------------------------------

    def addNeighbor(self,cell, nodes, style, domain, domains):

	geometrylist = list(nodes)

        faces = self.calculateFaces(style, geometrylist)
        self.domain.append(domain)
	self.neighbors.append(cell)
	self.geometry.append(geometrylist)
	self.faces.append(faces)
        celldata = [self.id, domain, cell, geometrylist, faces]
        self.celldata.append(celldata)

        """
        print self.id
        print self.position
        print self.domain
        print self.neighbors

        if (self.id==37):
            print self.id
            print self.geometry
            print self.faces
        """

#        celldata = [domain,cell,geometrylist,faces]

#        self.celldata.append(celldata)

# ---------------------------------------------

    def addCellData(self,cell, nodes, style, domain, domains):

	geometrylist = list(nodes)

        faces = self.calculateFaces(style, geometrylist)
        celldata = [domain,cell,geometrylist,faces]
        self.celldata.append(celldata)

#        print self.celldata
#        print "domain = " +str(self.domain) #+str(self.celldata.domain)
        
# ---------------------------------------------

    def printCellData(self):

        print "\n"
        print self.celldata.domain
        
###########################################        

    def addIndices(self,indices):
	self.indices=list(indices)

class CellData(object):

    def __init__(self, node):
        self.cell = node.neighbors
        self.domain = node.domain
        self.geometry = node.geometry
        self.faces = node.faces

class IndexedNodes(object):

    def __init__(self, ident, index, position,oldnodeid):
        self.id      = ident            # decimal number corresponding to indices
        self.indices = list(index)      # Logical grid location of the node
        self.position = list(position)  # Physical position of Node
	self.oldnodeid = oldnodeid

    def getID(self):
        return self.id

    def getOldNodeID(self):
        return self.oldnodeid

class Element(object):

    # Elements are the 'cells' of the Grid
    
    def __init__(self, ident, domain_index, nodes, neighbors):
        self.id = ident             # Unique Identity of element (global ID from Cubit)
        self.domain = domain_index  # The Domain this Element is found in
        self.nodes = list(nodes)    # Holds a list of node indices (not Global Indices)
        self.neighbors = list(neighbors)  # Holds a list of neighbor Element indices (not Global Indices)
        self.indices = []           # Holds the Global Indices of the element within the Grid
        self.located = False        # Tells us if this Element has been located within the grid

    # Returns the Global ID of the Element
    def getID(self):
        return self.id

    def getDomainIndex(self):
        return self.domain

    # Sets the Global Indices of the Element 
    # (Returns True if set, or False if Element already has a Global Index)
    def setIndices(self, globalIndices):
        if not self.located:
            self.indices = list(globalIndices)
            self.located = True
            return True
        else:
            return False

    # Return the Global Indices of the Element
    def getIndices(self):
        return list(self.indices)


class Domain(object):

    # Domains are groups of Elements within the Grid
    
    def __init__(self, domain_id, domain_name, primitive):
        self.id = domain_id             # Unique Global ID for Domain
        self.name = domain_name         # Unique Global Name for Domain (designated in Cubit)
        self.primitive = primitive      # Primitive type: Line, Rectangle or Hexahedral
        self.numDims = 0                # Dimensions within Domain (1D=[Line], 2D=[Line,Rectangle], 3D=[Line,Rectangle,Hexahedron])
        self.numElements = []           # Describes the number of Elements in each Dimension within the Domain
        self.numGhostCells = []         # Number of Ghost Cells at each end of a Dimension within Domain
        self.startIndices = []          # The starting Global Indices of the Domain (i_min, j_min, k_min)
        self.endIndices = []            # The ending Global Indices of the Domain (i_max, j_max, k_max)
        self.elements = []              # A list of the Elements within the Domain
        self.domaindata = 0             # This is used to hold information to export in hdf5
        self.blockfacelist = []         # Defines the faces of the domain.  
                                        # Sequence is (xpos,xneg,ypos,yneg,zpos,zneg).  
                                        # Each face defined with domain corner nodes.
                                        # Face definition for 
                                        # xpos: (node_x00,node_xy0,node_xyz,node_x0z)
                                        # xneg: (node_000,node_0y0,node_0yz,node_00z)
                                        # ypos: (node_0y0,node_xy0,node_xyz,node_0yz)
                                        # yneg: (node_000,node_x00,node_x0z,node_00z)
                                        # zpos: (node_00z,node_x0z,node_xyz,node_0yz)
                                        # zpos: (node_000,node_x00,node_xy0,node_0y0)
        self.n000 = []                  # contains nodeID corresponding to node_000 position on the block
        self.nx00 = []                  # contains nodeID corresponding to node_x00 position on the block
        self.n0y0 = []                  # contains nodeID corresponding to node_0y0 position on the block
        self.n00z = []                  # contains nodeID corresponding to node_00z position on the block
        self.nxy0 = []                  # contains nodeID corresponding to node_xy0 position on the block
        self.n0yz = []                  # contains nodeID corresponding to node_0yz position on the block
        self.nx0z = []                  # contains nodeID corresponding to node_x0z position on the block
        self.nxyz = []                  # contains nodeID corresponding to node_xyz position on the block
        self.geometry = []              # contains the following list of corner id's identifying the block 
                                        # [n000,n0y0,n0yz,n00z,nx00,nxy0,nxyz,nx0z]
        self.blockfaceneighborlist=[]       # a list of lists giving the neighboring domain for each face given in blockfacelist
        self.blockindexed=False         # boolean specification set to True after indexing     
        self.maxidx = []                # vector specifying the size of the indexed block in each direction
        self.indexedblock = []          # result of indexBlock operation
        self.indexedblockfaceslice = [] # Slice of data for a given face, following the ordering of blockfacelist, 
                                        # derived from indexedblock.
        self.indexedblockwgc = []       # indexedblock with ghost cells
        
    # Returns the Global ID of the Domain
    def getID(self):
        return self.id

    # Returns the Domain name
    def getName(self):
        return self.name

    # Returns the Primitive
    def getPrimitive(self):
        return self.primitive

    # Adds an Element to the Domain
    def addElement(self, elementIndex):
        self.elements.append(elementIndex)

    # Sets the Ghost Cells for the Domain
    def setGhostCells(self, ghostCells):
        self.numGhostCells = list(ghostCells)

    def getDomainFaceList(self):
        return self.blockfacelist

    def getDomainNeighborList(self):
        return self.blockneighborlist

    def addBlockFaceList(self,blockfacelist):
	self.blockfacelist=blockfacelist

    def addN000(self,n000id):
	self.n000=n000id

    def addNx00(self,nx00id):
	self.nx00=nx00id

    def addN0y0(self,n0y0id):
	self.n0y0=n0y0id

    def addN00z(self,n00zid):
	self.n00z=n00zid

    def addNxy0(self,nxy0id):
	self.nxy0=nxy0id

    def addN0yz(self,n0yzid):
	self.n0yz=n0yzid

    def addNx0z(self,nx0zid):
	self.nx0z=nx0zid

    def addNxyz(self,nxyzid):
	self.nxyz=nxyzid

    def printDomainCornerNodes(self,nodes):

        print 'n000 id: (x,y,z) = '+str(self.n000)+':\t('+str(round(nodes[self.n000].position[0],1))+','+str(round(nodes[self.n000].position[1],1))+','+str(round(nodes[self.n000].position[2],1))+')'
        print '         (i,j,k) = ('+str(nodes[self.n000].indices[0])+', '+str(nodes[self.n000].indices[1])+', '+str(nodes[self.n000].indices[2])+')'
        print 'n0y0 id: (x,y,z) = '+str(self.n0y0)+':\t('+str(round(nodes[self.n0y0].position[0],1))+','+str(round(nodes[self.n0y0].position[1],1))+','+str(round(nodes[self.n0y0].position[2],1))+')'
        print '         (i,j,k) = ('+str(nodes[self.n0y0].indices[0])+', '+str(nodes[self.n0y0].indices[1])+', '+str(nodes[self.n0y0].indices[2])+')'
        print 'n0yz id: (x,y,z) = '+str(self.n0yz)+':\t('+str(round(nodes[self.n0yz].position[0],1))+','+str(round(nodes[self.n0yz].position[1],1))+','+str(round(nodes[self.n0yz].position[2],1))+')'
        print '         (i,j,k) = ('+str(nodes[self.n0yz].indices[0])+', '+str(nodes[self.n0yz].indices[1])+', '+str(nodes[self.n0yz].indices[2])+')'
        print 'n00z id: (x,y,z) = '+str(self.n00z)+':\t('+str(round(nodes[self.n00z].position[0],1))+','+str(round(nodes[self.n00z].position[1],1))+','+str(round(nodes[self.n00z].position[2],1))+')'
        print '         (i,j,k) = ('+str(nodes[self.n00z].indices[0])+', '+str(nodes[self.n00z].indices[1])+', '+str(nodes[self.n00z].indices[2])+')'
        print 'nx00 id: (x,y,z) = '+str(self.nx00)+':\t('+str(round(nodes[self.nx00].position[0],1))+','+str(round(nodes[self.nx00].position[1],1))+','+str(round(nodes[self.nx00].position[2],1))+')'
        print '         (i,j,k) = ('+str(nodes[self.nx00].indices[0])+', '+str(nodes[self.nx00].indices[1])+', '+str(nodes[self.nx00].indices[2])+')'
        print 'nxy0 id: (x,y,z) = '+str(self.nxy0)+':\t('+str(round(nodes[self.nxy0].position[0],1))+','+str(round(nodes[self.nxy0].position[1],1))+','+str(round(nodes[self.nxy0].position[2],1))+')'
        print '         (i,j,k) = ('+str(nodes[self.nxy0].indices[0])+', '+str(nodes[self.nxy0].indices[1])+', '+str(nodes[self.nxy0].indices[2])+')'
        print 'nxyz id: (x,y,z) = '+str(self.nxyz)+':\t('+str(round(nodes[self.nxyz].position[0],1))+','+str(round(nodes[self.nxyz].position[1],1))+','+str(round(nodes[self.nxyz].position[2],1))+')'
        print '         (i,j,k) = ('+str(nodes[self.nxyz].indices[0])+', '+str(nodes[self.nxyz].indices[1])+', '+str(nodes[self.nxyz].indices[2])+')'
        print 'nx0z id: (x,y,z) = '+str(self.nx0z)+':\t('+str(round(nodes[self.nx0z].position[0],1))+','+str(round(nodes[self.nx0z].position[1],1))+','+str(round(nodes[self.nx0z].position[2],1))+')'
        print '         (i,j,k) = ('+str(nodes[self.nx0z].indices[0])+', '+str(nodes[self.nx0z].indices[1])+', '+str(nodes[self.nx0z].indices[2])+')'

        return
        
        
    def printDomainCornerNodeIndices(self,nodes):

        print 'n000 id: = '+str(self.n000)+':\t(i,j,k) = ('+str(nodes[self.n000].indices[0])+', '+str(nodes[self.n000].indices[1])+', '+str(nodes[self.n000].indices[2])+')'
        print 'n0y0 id: = '+str(self.n0y0)+':\t(i,j,k) = ('+str(nodes[self.n0y0].indices[0])+', '+str(nodes[self.n0y0].indices[1])+', '+str(nodes[self.n0y0].indices[2])+')'
        print 'n0yz id: = '+str(self.n0yz)+':\t(i,j,k) = ('+str(nodes[self.n0yz].indices[0])+', '+str(nodes[self.n0yz].indices[1])+', '+str(nodes[self.n0yz].indices[2])+')'
        print 'n00z id: = '+str(self.n00z)+':\t(i,j,k) = ('+str(nodes[self.n00z].indices[0])+', '+str(nodes[self.n00z].indices[1])+', '+str(nodes[self.n00z].indices[2])+')'
        print 'nx00 id: = '+str(self.nx00)+':\t(i,j,k) = ('+str(nodes[self.nx00].indices[0])+', '+str(nodes[self.nx00].indices[1])+', '+str(nodes[self.nx00].indices[2])+')'
        print 'nxy0 id: = '+str(self.nxy0)+':\t(i,j,k) = ('+str(nodes[self.nxy0].indices[0])+', '+str(nodes[self.nxy0].indices[1])+', '+str(nodes[self.nxy0].indices[2])+')'
        print 'nxyz id: = '+str(self.nxyz)+':\t(i,j,k) = ('+str(nodes[self.nxyz].indices[0])+', '+str(nodes[self.nxyz].indices[1])+', '+str(nodes[self.nxyz].indices[2])+')'
        print 'nx0z id: = '+str(self.nx0z)+':\t(i,j,k) = ('+str(nodes[self.nx0z].indices[0])+', '+str(nodes[self.nx0z].indices[1])+', '+str(nodes[self.nx0z].indices[2])+')'

        return
        
        
    def printDomainCornerNodeRoundedPositions(self,nodes):

        print 'n000 id: (x,y,z) = '+str(self.n000)+':\t('+str(round(nodes[self.n000].position[0],1))+','+str(round(nodes[self.n000].position[1],1))+','+str(round(nodes[self.n000].position[2],1))+')'
        print 'n0y0 id: (x,y,z) = '+str(self.n0y0)+':\t('+str(round(nodes[self.n0y0].position[0],1))+','+str(round(nodes[self.n0y0].position[1],1))+','+str(round(nodes[self.n0y0].position[2],1))+')'
        print 'n0yz id: (x,y,z) = '+str(self.n0yz)+':\t('+str(round(nodes[self.n0yz].position[0],1))+','+str(round(nodes[self.n0yz].position[1],1))+','+str(round(nodes[self.n0yz].position[2],1))+')'
        print 'n00z id: (x,y,z) = '+str(self.n00z)+':\t('+str(round(nodes[self.n00z].position[0],1))+','+str(round(nodes[self.n00z].position[1],1))+','+str(round(nodes[self.n00z].position[2],1))+')'
        print 'nx00 id: (x,y,z) = '+str(self.nx00)+':\t('+str(round(nodes[self.nx00].position[0],1))+','+str(round(nodes[self.nx00].position[1],1))+','+str(round(nodes[self.nx00].position[2],1))+')'
        print 'nxy0 id: (x,y,z) = '+str(self.nxy0)+':\t('+str(round(nodes[self.nxy0].position[0],1))+','+str(round(nodes[self.nxy0].position[1],1))+','+str(round(nodes[self.nxy0].position[2],1))+')'
        print 'nxyz id: (x,y,z) = '+str(self.nxyz)+':\t('+str(round(nodes[self.nxyz].position[0],1))+','+str(round(nodes[self.nxyz].position[1],1))+','+str(round(nodes[self.nxyz].position[2],1))+')'
        print 'nx0z id: (x,y,z) = '+str(self.nx0z)+':\t('+str(round(nodes[self.nx0z].position[0],1))+','+str(round(nodes[self.nx0z].position[1],1))+','+str(round(nodes[self.nx0z].position[2],1))+')'

        return
        
    def printDomainCornerNodeIDs(self,nodes):
        
        print 'domainID = '+str(self.id)

        print 'n000 id: '+str(self.n000)
        print 'n0y0 id: '+str(self.n0y0)
        print 'n0yz id: '+str(self.n0yz)
        print 'n00z id: '+str(self.n00z)
        print 'nx00 id: '+str(self.nx00)
        print 'nxy0 id: '+str(self.nxy0)
        print 'nxyz id: '+str(self.nxyz)
        print 'nx0z id: '+str(self.nx0z)
        return
        
    def printDomainCornerGeometryNodeIDs(self,nodes):
        
        print 'domainID: '+str(self.id)+' - geometry: '+str(self.geometry)
        return


    def findDomainNeighbor(self, nodes):

        domainfacelist         = self.blockfacelist
        domainfaceneighborlist = self.blockfaceneighborlist

#        print 'findDomainNeighbor:          domainfacelist = ' + str(domainfacelist)
#        print 'findDomainNeighbor:  domainfaceneighborlist = ' + str(domainfaceneighborlist)

        faces = range(0,6)

        print 'findDomainNeighbor:                domainID = '+str(self.id)
        for face in faces:
            print 'findDomainNeighbor:                    face = '+str(face)
            currentface = domainfacelist[face]
            
            setlist = []

            for node in range(len(currentface)):
                nodeID = currentface[node]
                domainlist = deepcopy(nodes[nodeID].domain)
                domainlist.remove(self.id)
                setlist.append(set(domainlist))
 
            tempfacedomainlist = list(set.intersection(*setlist))

                                    
            self.blockfaceneighborlist.append(tempfacedomainlist)
            domainfaceneighborlist = self.blockfaceneighborlist
            print 'findDomainNeighbor:  domainfaceneighborlist = ' + str(domainfaceneighborlist)
            print 'findDomainNeighbor:  pause'
            
        print 'findDomainNeighbor:  domain ' +str(self.id) + ' facelist = ' +str(domainfacelist)
        print 'findDomainNeighbor:  domainfaceneighborlist = ' + str(domainfaceneighborlist)

        return

    def generateBlockCornerNodes(self, indexednodes, maxidx, style):

        cornerNodeIDs = getCornerNodeIDs(indexednodes,maxidx)

        node_000 = cornerNodeIDs[0]
        node_x00 = cornerNodeIDs[1]
        node_0y0 = cornerNodeIDs[2]
        node_xy0 = cornerNodeIDs[3]
        node_00z = cornerNodeIDs[4]
        node_x0z = cornerNodeIDs[5]
        node_0yz = cornerNodeIDs[6]
        node_xyz = cornerNodeIDs[7]

        self.geometry = [node_000,node_0y0,node_0yz,node_00z,node_x00,node_xy0,node_xyz,node_x0z]

        print 'generateBlockFaces:  geometry = ' +str(self.geometry)

        self.addN000(node_000)
        self.addNx00(node_x00)
        self.addN0y0(node_0y0)
        self.addN00z(node_00z)
        self.addNxy0(node_xy0)
        self.addN0yz(node_0yz)
        self.addNx0z(node_x0z)
        self.addNxyz(node_xyz)

        return


    def generateBlockFaces(self, indexednodes, maxidx, style):

        self.generateBlockCornerNodes(indexednodes, maxidx, style)

        node_000 = self.n000
        node_x00 = self.nx00
        node_0y0 = self.n0y0
        node_00z = self.n00z
        node_xy0 = self.nxy0
        node_0yz = self.n0yz
        node_x0z = self.nx0z
        node_xyz = self.nxyz

        face_xpos = [node_x00,node_x0z,node_xyz,node_xy0]
        face_xneg = [node_000,node_00z,node_0yz,node_0y0]
        
        face_ypos = [node_0y0,node_0yz,node_xyz,node_xy0]
        face_yneg = [node_000,node_00z,node_x0z,node_x00]
        
        face_zpos = [node_00z,node_0yz,node_xyz,node_x0z]
        face_zneg = [node_000,node_0y0,node_xy0,node_x00]

        blockfacelist = [face_xpos, face_xneg, face_ypos, face_yneg, face_zpos, face_zneg]
        
        print 'domain '+str(self.id)+' facelist = '+str(blockfacelist)

        self.addBlockFaceList(blockfacelist)

        return 

    def calculateFacesFromGeometry(self, style, geometrylist):
        
        """
        This function calculates the face nodelists given a geometry
        """

        if (style == 'Hexahedron'):
            face1pindices,face1nindices,face2pindices,face2nindices,face3pindices,face3nindices = cellHexahedronFaceIndices()
            face1p = [geometrylist[i] for i in face1pindices]
            face1n = [geometrylist[i] for i in face1nindices]
            face2p = [geometrylist[i] for i in face2pindices]
            face2n = [geometrylist[i] for i in face2nindices]
            face3p = [geometrylist[i] for i in face3pindices]
            face3n = [geometrylist[i] for i in face3nindices]
            faces = (face1p,face1n,face2p,face2n,face3p,face3n)
        elif (style == 'Rectangle'):
            face1nindices = (0,1)
            face2nindices = (0,3)
            face1pindices = (3,2)
            face2pindices = (1,2)
            face1p = [geometrylist[i] for i in face1pindices]
            face1n = [geometrylist[i] for i in face1nindices]
            face2p = [geometrylist[i] for i in face2pindices]
            face2n = [geometrylist[i] for i in face2nindices]
            faces = (face1p,face1n,face2p,face2n)

        self.addBlockFaceList(faces)
    
        return


    def printDomainCornerNodePositions(self,nodes):
    
        print 'n000 position, indices = '+str(nodes[self.n000].position)
        print 'nx00 position, indices = '+str(nodes[self.nx00].position)
        print 'n0y0 position, indices = '+str(nodes[self.n0y0].position)
        print 'n00z position, indices = '+str(nodes[self.n00z].position)
        print 'nxy0 position, indices = '+str(nodes[self.nxy0].position)
        print 'n0yz position, indices = '+str(nodes[self.n0yz].position)
        print 'nx0z position, indices = '+str(nodes[self.nx0z].position)
        print 'nxyz position, indices = '+str(nodes[self.nxyz].position)

        return
            
    def getAndStoreIndexedBlockFaceSlices(self):
        
        numpyarray = self.indexedBlock
        ghostcells = self.numGhostCells
                
        for facepos in range(6):
            
            ghostcellslicethickness = ghostcells[facepos]       
            localghostcells = range(ghostcellslicethickness)

            if ((facepos==0)or(facepos==1)):            
                direction = 0
            elif ((facepos==2)or(facepos==3)):            
                direction = 1
            elif ((facepos==4)or(facepos==5)):            
                direction = 2
                
            print 'self.maxidx[direction] = '+str(self.maxidx[direction])
            print 'localghostcells = '+str(localghostcells)
            
            if ((facepos==0)or(facepos==2)or(facepos==4)):
                ghostcellrange = [self.maxidx[direction]-lgc for lgc in localghostcells]
            else:
                ghostcellrange = localghostcells

            if (direction==0):
                ii,jj,kk = ghostcellrange,s_[:],s_[:]
            elif (direction==1):
                ii,jj,kk = s_[:],ghostcellrange,s_[:]
            elif (direction==2):
                ii,jj,kk = s_[:],s_[:],ghostcellrange
            
            ghostcellslice = numpyarray[ii][jj][kk]
            self.indexedblockfaceslice.append(ghostcellslice)
        
        return
        
    def constructNodeArrayWithGhostCells(self,slices):
        
        numpynodearraycore = copy(self.indexedBlock)
        maxidx = self.maxidx
        ghostcells = self.numGhostCells # vector of the number of ghostcells at the end of each direction
        numofparams = 5
        
        maxidxwgc[0] = maxidx[0]+ghostcells[0]+ghostcells[1]        
        maxidxwgc[1] = maxidx[1]+ghostcells[2]+ghostcells[3]        
        maxidxwgc[2] = maxidx[2]+ghostcells[4]+ghostcells[5]        
        
        self.maxidxwgc = maxidxwgc
        
        ghostcornerblock = ones((ghostcells[0],ghostcells[2],ghostcells[4],numofparams))
        ghostedgeblockx = ones((maxidx[0],ghostcells[2],ghostcells[4],numofparams))
        ghostedgeblocky = ones((ghostcells[0],maxidx[1],ghostcells[4],numofparams))
        ghostedgeblockz = ones((ghostcells[0],ghostcells[2],maxidx[2],numofparams))

        ghostslicexpos = slices[0]
        ghostslicexneg = slices[1]
        ghostsliceypos = slices[2]
        ghostsliceyneg = slices[3]
        ghostslicezpos = slices[4]
        ghostslicezneg = slices[5]

        edgeblockx = concatenate((ghostcornerblock,ghostedgeblockx,ghostcornerblock),axis=0) 
        edgeblocky = concatenate((ghostcornerblock,ghostedgeblocky,ghostcornerblock),axis=1) 
        edgeblockz = concatenate((ghostcornerblock,ghostedgeblockz,ghostcornerblock),axis=2) 
        
        block1 = concatenate((ghostslicexneg,numpynodearraycore, ghostslicexpos),axis=0)
        block2 = concatenate((edgeblockz,ghostsliceypos,edgeblockz),axis=0)
        block3 = concatenate((edgeblockz,ghostsliceyneg,edgeblockz),axis=0)

        block4 = concatenate((block3,block1,block2),axis=1)
        
        block5 = concatenate((edgeblocky,ghostslicezpos,edgeblocky),axis=0)
        block6 = concatenate((edgeblocky,ghostslicezneg,edgeblocky),axis=0)
        
        block7 = concatenate((ghostcornerblock,edgeblockx,ghostcornerblock),axis=0)
        
        block8 = concatenate((block7,block5,block7),axis=1)       
        block9 = concatenate((block7,block6,block7),axis=1)     
        
        finalblock = concatenate((block9, block4, block8), axis=2)
        
        self.blockindexedwgc = finalblock
        
        return
        
    def reindexEnlargedNumpyNodeArray(self,nodes):

        enlargednumpynodearray= self.blockindexedwgc

        imax = self.maxidxwgc[0]+1
        jmax = self.maxidxwgc[1]+1
        kmax = self.maxidxwgc[2]+1

        for i in range(imax):
            for j in range(jmax):
                for k in range(kmax):

                    newindexednodeID = kmax*jmax*i + kmax*j + k
                    enlargednumpynodearray[i][j][k].id = newindexednodeID

                    enlargednumpynodearray[i][j][k].indices[0] = i
                    enlargednumpynodearray[i][j][k].indices[1] = j
                    enlargednumpynodearray[i][j][k].indices[2] = k
                    
                    nodeid_old = enlargednumpynodearray[i][j][k].oldnodeid
                    
                    nodes[nodeid_old].indices[0] = i
                    nodes[nodeid_old].indices[1] = j
                    nodes[nodeid_old].indices[2] = k
  
        self.blockindexedwgc = enlargednumpynodearray

        return
        
    def translateEnlargedNumpyNodeArrayBackToPythonNodeArray(self):

        enlargednumpynodearray= self.blockindexedwgc

        imax = self.maxidxwgc[0]+1
        jmax = self.maxidxwgc[1]+1
        kmax = self.maxidxwgc[2]+1

        for i in range(imax):
            for j in range(jmax):
                for k in range(kmax):

                    indices          = [i,j,k]
                    newindexednodeID = kmax*jmax*i + kmax*j + k

                    oldnodeid = enlargednumpynodearray[i][j][k].oldnodeid
                    position  = enlargednumpynodearray[i][j][k].position
                    
                    nodearray.append(IndexedNodes(newindexednodeID,indices,position,oldnodeid))

        sortID(nodearray)
    
        self.finalnodearray = nodearray

        return   
        
    def finalBlockConstruction(self,nodes,domains):

        slices = []
    
        for face in range(6):

            print 'finalBlockConstruction: current face is '+str(face)

            nextdomainlistlength = len(self.blockfaceneighborlist[face])
            facepresentflag      = (nextdomainlistlength!=0)

            if (facepresentflag):
                nextdomainIDlist  = self.blockfaceneighborlist[face]
                nextdomainID      = nextdomainIDlist[0]
                nextdomain        = domains[nextdomainID]
                if (face == 0):
                    nextdomainface=1
                elif (face == 1):
                    nextdomainface=0
                elif (face == 2):
                    nextdomainface=3
                elif (face == 3):
                    nextdomainface=2
                elif (face == 4):
                    nextdomainface=5
                elif (face == 5):
                    nextdomainface=4
                                
                slices.append(nextdomain.indexedblockfaceslice[nextdomainface])
            
            else:
                slices.append(self.indexedblockfaceslice[face])
  
        self.constructNodeArrayWithGhostCells(slices)
        
        self.reindexEnlargedNumpyNodeArray(nodes)
        
        self.translateEnlargedNumpyNodeArrayBackToPythonNodeArray()
        
        
           
        domainID = self.id
        print "finalBlockConstruction: Saving Block..." +str(domainID)        
        outfilename = filename + "_block_" +str(domainID)+ ".inp"    
        writeASCIIFile(outfilename, self.finalnodearray, nodes, style)
        print "finalBlockConstruction: Written to " + outfilename

        return          
                
def searchForNextNode(domains, oldnocells, nodes, nnID, nodeID, cellID, domainID, idv, nfv, nnv, style, celllist, oddslice,firstblockpath):

    celldelta, sharedcells     = cellDelta(oldnocells,nodes,nnID,domainID)
    nodevisited                = 0
    previousDomainID           = -1
    
#    celllist                   = nodes[nodeID].findSharedCellsForSingleDomain(domainID)
    
    if (cellID!=-1):
        oldcellID              = cellID
        
    oldnodecellposition,cellID = findCell(nnID, oldcellID, nfv, nodes, idv, style, celldelta, nodevisited, domainID,firstblockpath) 
#    fpv, nIDp, cellIDposition  =               defineFaceandCellPosition(domains, nodes, nnID, cellID,           domainID,previousDomainID,firstblockpath,style)
    fpv, nIDp, cellIDposition  = defineFaceandCellPosition(domains, nodes, nnID, cellID, domainID,previousDomainID,firstblockpath,style)
#    fpv, nIDp, cellIDposition  = defineFaceandCellPosition(domains, nodes,nodeID,cellID,domainID,previousDomainID,firstblockflag,style)
#    print 'searchForNextNode:  fpv = '+str(fpv)+' and idv = '+str(idv)
    nfv                        = buildNfv(fpv,style)
    nnv                        = buildNnv(nfv,nIDp,nodes,nnID,cellID,cellIDposition,domainID)
    oldnodeID                  = nnID
    nnID                       = getNextNode(nnv,idv)
    
#    veryoldnocells             = oldnocells
#    oldcelldelta               = celldelta
    oldnocells                 = sharedcells 
    
    nodeID                     = nnID

    return nodeID, oldnocells, nnID, oldnodeID, nnv, nfv, fpv, nIDp, cellID,nodevisited,celldelta

def findDomainCornerNode(domains, nodes, nodeID, cellID, style, domainID, idv, firstblockpath):

    nfv, nnv = setNfvandNnv(domains, nodes,nodeID,cellID,style,domainID,firstblockpath)

    nnID = getNextNode(nnv,idv)

#    firstnodeflag = False
#    firstendofrow = True
#    lengthofi  = 0
    oldnocells = 1
    oldnodeID  = nodeID
    oddslice = True

    celllist  = nodes[nodeID].findSharedCellsForSingleDomain(domainID)
    nodeID = nnID

    numberofnodestoscan = range(len(nodes)-1)

    firstblockpath = True

    for node in numberofnodestoscan:

        nodeID, oldnocells, nnID, oldnodeID, nnv, nfv, fpv, nIDp, cellID,nodevisited,celldelta = searchForNextNode(domains, oldnocells, nodes, nnID, nodeID, cellID, domainID, idv, nfv, nnv, style, celllist,oddslice,firstblockpath)

        thisnode = nodes[nodeID]

#        singledomainneighborlist = thisnode.findSharedCellsForSingleDomain(domainID)
        numofneighbors= thisnode.findNumofSharedCells(domainID)

        if (numofneighbors==1):
            cornernodeID = nodeID
#            print 'cornernodeID = '+str(cornernodeID)
            break

#    print 'returning from findDomainCornerNode'

    return cornernodeID, cellID

def buildOppositeCellFace(nodes, nodeID, cellID,sharedfacenodelist,style):

        oppositeface = [0, 0, 0, 0]
	idv = [0, 0, 0]

#        geometry = nodes[nodeID].geometry[cellIDposition]

#        print 'buildOppositeCellFace:  blocksharedfacenodelist = ' +str(blocksharedfacenodelist)

        for nodeposition in range(len(blocksharedfacenodelist)):

#            nodeID = sharedfacenodelist[nodeposition]

            idv[0] = 1
            idv[1] = 0
            idv[2] = 0

#            node             = nodes[nodeID]
#            nodeneighborlist = node.findSharedCellsForSingleDomain(domainID)
 #           oldcellID        = nodeneighborlist[0]

            for dirck in range(len(idv)):
                     
#                cellcornernodeID = geometry[geometrylistposition]

                for nodeIDcheck in blocksharedfacenodelist:

                    nodematch = (nodeIDcheck==newcornernodeID)

                    if (nodematch):
                        if (idv[0]==1):
                            idv[0]=0
                            idv[1]=1
                        elif (idv[1]==1):
                            idv[1]=0
                            idv[2]=1
                        elif (idv[2]==1):
                            idv[2]=0
                            idv[0]=1
                        break


            if not(nodematch):
                oppositeface[nodeposition]=newcornernodeID
                

#        print 'oppositeface = '+str(oppositeface)

        return oppositeface        


    
def buildOppositeFace(domains, blocksharedfacenodelist,domainID,nodes,style, firstblockpath):
    
        """
        Overview:
            
        To build the opposite face for the block, the basic idea is to take the 
        blocksharedfacenodelist and find the matching nodes for each corner for the 
        requested domain.  
        
        Detail:
            
        Each node from the blocksharedfacenodelist is taken as the starting point.
        then each direction is specified and the nodes object is scanned via 
        findDomainCornerNode to find the corner that moving through the block in
        the specified direction.  A newcornernodeID and cell ID are returned.
        
        As yet it is unknown which of the 3 directions is actually being specified, 
        i.e. x, y, or z.  However, it is not necessary to know at this point what the
        direction actually is.  All that is required is finding that cornernode on the 
        opposite face.
        
        Since the direction is unknown, the newcornernodeID is compared against the 
        blocksharedfacenodelist.  This is because the odds are 2 out of 3 that the 
        newcornernodeID determined is actually on the blocksharedfacenodelist, and 
        not on the opposite face.  If the test is positive, and the newcornernodeID 
        is in fact on the blocksharedfacenodelist, then the direction is incremented 
        and the test is repeated. If, however, the test happens to fail, then the 
        assignment is made of the newcornernodeID to the oppositeface (nodelist). 
        
        The cost of this approach is that it requires scanning the edge of the block
        a total of 3 times to find one node on the opposite face.

        """
    
        oppositeface = [0, 0, 0, 0]
        idv = [0, 0, 0]
    
#        print 'buildOppositeFace:  blocksharedfacenodelist = ' +str(blocksharedfacenodelist)

        for nodeposition in range(len(blocksharedfacenodelist)):

            nodeID = blocksharedfacenodelist[nodeposition]

            idv[0] = 1
            idv[1] = 0
            idv[2] = 0

            node             = nodes[nodeID]
            nodeneighborlist = node.findSharedCellsForSingleDomain(domainID)
            oldcellID        = nodeneighborlist[0]

            for dirck in range(len(idv)):
                     
                newcornernodeID, cellID = findDomainCornerNode(domains, nodes,nodeID,oldcellID,style,domainID, idv, firstblockpath) 

                for nodeIDcheck in blocksharedfacenodelist:

                    nodematch = (nodeIDcheck==newcornernodeID)

                    if (nodematch):
                        if (idv[0]==1):
                            idv[0]=0
                            idv[1]=1
                        elif (idv[1]==1):
                            idv[1]=0
                            idv[2]=1
                        elif (idv[2]==1):
                            idv[2]=0
                            idv[0]=1
                        break


            if not(nodematch):
                oppositeface[nodeposition]=newcornernodeID
                

#        print 'oppositeface = '+str(oppositeface)

        return oppositeface        

#"""
def calculateDomainCorners(style, blocksharedfacenodelist, blockoppositefacenodelist, blocksharedfaceposition,domains,nextDomainID,currentDomainID):
    
    """
    Overview:
        
    The purpose of this module is to find the corners of the nextDomainID, and the 
    corresponding geometry listing of these cornernodes from the 
    blocksharedfacenodelist and blockoppositefacenodelist.  The geometry vector for
    the next block (specified with nextDeomainID) is reconstructed.  Once done, the 
    corner nodes are stored, and the geometry vector is returned.
    """

    cD = domains[currentDomainID]

    cD_n000id = cD.n000
    cD_nx00id = cD.nx00
    cD_n0y0id = cD.n0y0
    cD_n00zid = cD.n00z
    cD_nxy0id = cD.nxy0
    cD_n0yzid = cD.n0yz
    cD_nx0zid = cD.nx0z
    cD_nxyzid = cD.nxyz

    cD_geometry = [cD_n000id, cD_n0y0id, cD_n0yzid, cD_n00zid, cD_nx00id, cD_nxy0id, cD_nxyzid, cD_nx0zid]
    
    knownface = cD.blockfacelist[blocksharedfaceposition]

    face1pindices,face1nindices,face2pindices,face2nindices,face3pindices,face3nindices = cellHexahedronFaceIndices()

    indices = [face1pindices,face1nindices,face2pindices,face2nindices,face3pindices,face3nindices]

    facenodelist = [cD_geometry[i] for i in indices[blocksharedfaceposition]]    

    matchednodes = 0
    
    bsfnodepos_indices = []
    
    for facenodepos in range(len(facenodelist)):
        facenodeid = facenodelist[facenodepos]
        for bsfnodepos in range(len(blocksharedfacenodelist)):
            bsfnodeid = facenodelist[bsfnodepos]
            if (bsfnodeid == facenodeid):
                matchednodes = matchednodes+1
                bsfnodepos_indices.append(bsfnodepos)
                
    bsfnodepos_indices = [0, 1, 2, 3]
    print 'bsfnodepos_indices = ' +str(bsfnodepos_indices)
    
                
    if (matchednodes!=len(facenodelist)):
        print 'calculateDomainFacesError:  matchednodes fails test.'
    else:
        if (blocksharedfaceposition==0):
            nextDomain_blocksharedfaceposition = 1
        elif (blocksharedfaceposition==2):
            nextDomain_blocksharedfaceposition = 3
        elif (blocksharedfaceposition==4):
            nextDomain_blocksharedfaceposition = 5
        elif (blocksharedfaceposition==1):
            nextDomain_blocksharedfaceposition = 0
        elif (blocksharedfaceposition==3):
            nextDomain_blocksharedfaceposition = 2
        elif (blocksharedfaceposition==5):
            nextDomain_blocksharedfaceposition = 4
        
    nextDomain_blockoppositefaceposition = blocksharedfaceposition

    nextDomain_bsf_geometryindices = indices[nextDomain_blocksharedfaceposition]
    nextDomain_bof_geometryindices = indices[nextDomain_blockoppositefaceposition]
    
    bofnodepos_indices = bsfnodepos_indices

    nD_sharedface   = [blocksharedfacenodelist[i] for i in bsfnodepos_indices]
    
#    nD_sharedface = [0,0,0,0]
#    for ipos in range(len(bsfnodepos_indices)):
#        i = bsfnodepos_indices[ipos]
#        nD_sharedface[ipos]=blocksharedfacenodelist[i]
    
    nD_oppositeface = [blockoppositefacenodelist[i] for i in bofnodepos_indices]
    
    print 'nD_sharedface = '+str(nD_sharedface)
    print 'nD_oppositeface = '+str(nD_oppositeface)
    
    nD_geometry = [0,0,0,0,0,0,0,0]

    for i in range(len(nextDomain_bsf_geometryindices)):
        nextDomain_bsf_geometryindex = nextDomain_bsf_geometryindices[i]
        nextDomain_bof_geometryindex = nextDomain_bof_geometryindices[i]
        nD_geometry[nextDomain_bsf_geometryindex] = nD_sharedface[i]
        nD_geometry[nextDomain_bof_geometryindex] = nD_oppositeface[i]
#        print 'nD_sharedface['+str(i)+'] = '+str(nD_sharedface[i])
#        print 'nD_geometry['+str(nextDomain_bsf_geometryindex)+'] = '+str(nD_geometry[nextDomain_bsf_geometryindex])

    print 'nD_geometry = ' +str(nD_geometry)
    
    if (style == 'Hexahedron'):
        node_000 = nD_geometry[0]
        node_0y0 = nD_geometry[1]
        node_0yz = nD_geometry[2]
        node_00z = nD_geometry[3]
        node_x00 = nD_geometry[4]
        node_xy0 = nD_geometry[5]
        node_xyz = nD_geometry[6]
        node_x0z = nD_geometry[7]
        domains[nextDomainID].addN000(node_000)
        domains[nextDomainID].addNx00(node_x00)
        domains[nextDomainID].addN0y0(node_0y0)
        domains[nextDomainID].addN00z(node_00z)
        domains[nextDomainID].addNxy0(node_xy0)
        domains[nextDomainID].addN0yz(node_0yz)
        domains[nextDomainID].addNx0z(node_x0z)
        domains[nextDomainID].addNxyz(node_xyz)
    else:
        print 'hexahedron not selected, and is for now the only option'
#   
#    cD_geometry = [cD_n000id, cD_n0y0id, cD_n0yzid, cD_n00zid, cD_nx00id, cD_nxy0id, cD_nxyzid, cD_nx0zid]
#

    return nD_geometry
#"""



class Face(object):

    # The Face class just defines an Elements facial Nodes for building the Connectivity of the Grid
    
    def __init__(self, n1, n1_f, n2, n2_f, nodes, n1dom, n2dom):
        self.n1 = n1            # Element index on side '1' of face
        self.n1_f = n1_f        # Face ID of Face for Element 1
        self.n2 = n2            # Element index on side '2' of face
        self.n2_f = n2_f        # Face ID of Face for Element 2
        self.n1_d = n1dom       # Domain ID for Element 1
        self.n2_d = n2dom       # Domain ID for Element 2
        self.nodes = nodes      # Nodes found on Face
        self.nodes.sort()

    # Compares Face 'self' with face 'other' (Returns '-1' if 'self' is smaller, 1 if 'self' is larger, 0 if equal) 
    def compareFaces(self, other):
        for i in range(len(self.nodes)):
            if self.nodes[i] > other.nodes[i]:
                return 1
            elif self.nodes[i] < other.nodes[i]:
                return -1
        return 0
    
########################################################################################
#### File IO

# Sorts a list by the item's IDs
def sortID(itemlist):

    # I grabbed this from somwhere... its fast and useful
    # I'm not sure how fast it is compared to a merge sort in this case
    
    def quick_sort(list2):
        quick_sort_r(list2, 0, len(list2) - 1)
 
    # quick_sort_r, recursive (used by quick_sort)
    def quick_sort_r(list2 , first, last):
        if last > first:
            pivot = partition(list2, first, last)
            quick_sort_r(list2, first, pivot - 1)
            quick_sort_r(list2, pivot + 1, last)

    def partition(list2, first, last):
        sred = (first + last)/2
        if list2[first].getID() > list2[sred].getID():
            list2[first], list2[sred] = list2[sred], list2[first]  # swap
        if list2[first].getID() > list2[last].getID():
            list2[first], list2[last] = list2[last], list2[first]  # swap
        if list2[sred].getID() > list2[last].getID():
            list2[sred], list2[last] = list2[last], list2[sred]    # swap
            
        list2 [sred], list2 [first] = list2[first], list2[sred]    # swap
        pivot = first
        i = first + 1
        j = last
     
        while True:
            while i <= last and list2[i].getID() <= list2[pivot].getID():
                i += 1
            while j >= first and list2[j].getID() > list2[pivot].getID():
                j -= 1
            if i >= j:
                break
            else:
                list2[i], list2[j] = list2[j], list2[i]  # swap
        list2[j], list2[pivot] = list2[pivot], list2[j]  # swap
        return j
            
    quick_sort(itemlist)

# Runs a binary search through a sorted list and returns the index of the item
def binaryIDSearch(itemlist, itemid):
    start = 0;
    end = len(itemlist)
    while start < end:
        mid = (end + start) / 2
        midid = itemlist[mid].getID()
        if midid < itemid:
            start = mid
        elif midid > itemid:
            end = mid
        else:
            return mid
    return -1

# Returns an item's index in a list given an Unique ID
def getIndexFromID(itemlist, itemid):
    itemindex = itemid - 1
    if itemindex < len(itemlist):
        if itemlist[itemindex].getID() == itemid:
            return itemindex
    return binaryIDSearch(itemlist, itemid)

# Reads an 'Abacus' File and returns the Nodes, Elements and Domains of the Grid
def readFile(filename, ghostCells):

    # Open the file, grab its contents and close the file
    f = open(filename, "r")
    lines = f.readlines()
    f.close()

    # Now we run through the lines of the file, pulling out the information we need

    # First we set some constants:
    style = ""              # This defines the polygon primative: ""=Error, "Lines", "Rectangles", "Hexahedrons"
    elements = []           # This is a list of the Elements found in the File
    nodes = []              # This is a list of Nodes found in the File
    domains = []            # This is a list of Domains in the Grid
    n_faces = 0             # This is a placeholder for defining how many Faces are needed for each Element
    spacing = "\t"          # This is for printing to screen
    numlines = len(lines)   # This is the number of Lines read from the File
    i = 0                   # This is the index of the current Line
    while True:
        # There are two keys we are looking for in the Lines: Nodes and Elements

        # Nodes:
        if(lines[i] == "********************************** N O D E S **********************************\n"):
            # Skip a Line
            i += 2
            
            # All nodes are in one big group
            while True:

                # Check if this Line contains useful information, else we end our search for Nodes
                if lines[i][0] == "*" or lines[i] == "":
                    break
                
                # First we Evaluate the Line to build a tuple
                postuple = eval(lines[i])
                
                # Then we create a Node from the tuple
                nodes.append(Node(postuple[0], postuple[1:]))
                i += 1

        # Elements:
        elif(lines[i] == "********************************** E L E M E N T S ****************************\n"):
            
            i += 1
            domainIndex = -1   # This will be the ID of each domain through reconstruction
            elementIndex = 0 # This is used to locate the Elements based on their Index (muuuuuch faster)
            
            # We run through the Lines to collect the Domains
            while True:
                # Check if this Line contains useful data, else we end our search for Elements
                if lines[i][:2] == "**" or lines[i] == "":
                    break
                
                # If this line designates a new Domain of Elements:
                if lines[i][:2] == "*E":
                    # Increment the domainIndex
                    domainIndex += 1
                    # Check what kind of Domain this is
                    if lines[i][15:20] == "C3D8R":      # Hexahedron
                        n_faces = 6
                        domain_name = lines[i][28:].splitlines()[0]
                        style = "Hexahedron"
                    elif lines[i][15:18] == "S4R":      # Rectangles
                        n_faces = 4
                        domain_name = lines[i][26:].splitlines()[0]
                        style = "Rectangle"
                    elif lines[i][15:18] == "B21":      # Lines
                        n_faces = 2
                        domain_name = lines[i][26:].splitlines()[0]
                        style = "Line"
                    else:
                        print "This grid converter does not work with the Element primitive in the input file"
                        return [],[],[],0

                    # Now we append this new Domain to our list of Domains
                    domains.append(Domain(domainIndex, domain_name, style))
                    domains[-1].setGhostCells(ghostCells)
#                    domainIndex += 1
                    
                else:
                    # This line designates a new Element within a Domain

                    # First we grab the information from the Line
                    elementinfo = eval(lines[i])

                    # Next we append a new Element to the list of Elements
                    elements.append(Element(elementinfo[0], domainIndex, elementinfo[1:], [0 for j in range(n_faces)]))

                    # We need to add this Element's index to the Domain it belongs to
                    domains[-1].addElement(elementIndex)

                    # Now increment the Element Index
                    elementIndex += 1
                
                i += 1
        
        else:
            # Increment to the next Line if nothing has happened
            i += 1

        # This is our exit clause for when we run out of file
        if(i >= numlines):
            break

    # Print out what we found:
    print spacing + "Nodes Found: " + str(len(nodes))
    print spacing + "Domains Found: " + str(len(domains))
    print spacing + "Elements Found: " + str(len(elements))

    # The last thing to do before we return this data is to
    # sort the Nodes by their ID and map the Elements to
    # the Node Indices instead of the Node IDs

    print spacing + "Sorting Nodes..."
    sortID(nodes)

    print spacing + "Nodes Sorted"
    print spacing + "Reindexing Element Nodes..."
    for i in range(len(elements)):
        nodeset = elements[i].nodes
        for j in range(len(nodeset)):
            nodeIndex = getIndexFromID(nodes,nodeset[j])
            elements[i].nodes[j] = nodeIndex
        #print "Rebuilding Element " + str(i) + ": " + str(elements[i].nodes)
    print spacing + "Element Nodes Reindexed"
    
    for node in range(len(nodes)):
	nodes[node].id=node
    
    return nodes, elements, domains, style

# The following functions take a list of elements and builds a list of faces from them (each face is only connected to a single Element initially)

def buildLineFaces(nodes,currentnode):
    faces = []
    face_id = [[0],[1]] # This is untested, but probably correct
    for i in range(len(elements)):
        faces.append(Face(i,0,-1,0,[elements[i].nodes[face_id[0][0]]],elements[i].getDomainIndex(), -1))
        faces.append(Face(i,1,-1,0,[elements[i].nodes[face_id[1][0]],elements[i].getDomainIndex(), -1]))
    return faces

def buildRectangleFaces(nodes,currentnode):
    faces = []
    face_id = [[0,1],[1,2],[2,3],[3,0]] # This is untested but probably correct
    for i in range(len(nodes[currentnode].neighbor)):
        
        faces.append(Face(i,0,-1,0,[elements[i].nodes[face_id[0][0]],elements[i].nodes[face_id[0][1]]],elements[i].getDomainIndex(), -1))
        faces.append(Face(i,1,-1,0,[elements[i].nodes[face_id[1][0]],elements[i].nodes[face_id[1][1]]],elements[i].getDomainIndex(), -1))
        faces.append(Face(i,2,-1,0,[elements[i].nodes[face_id[2][0]],elements[i].nodes[face_id[2][1]]],elements[i].getDomainIndex(), -1))
        faces.append(Face(i,3,-1,0,[elements[i].nodes[face_id[3][0]],elements[i].nodes[face_id[3][1]]],elements[i].getDomainIndex(), -1))
    return faces

def buildHexahedronFaces(nodes,currentnode):
    faces = []
    face_id = [[0,1,2,3],[4,5,6,7],[0,1,4,5],[2,3,6,7],[1,2,5,6],[0,3,4,7]]
    for i in range(len(elements)):
        faces.append(Face(i,0,-1,0,[elements[i].nodes[face_id[0][0]],elements[i].nodes[face_id[0][1]],elements[i].nodes[face_id[0][2]],elements[i].nodes[face_id[0][3]]],elements[i].getDomainIndex(), -1))
        faces.append(Face(i,1,-1,0,[elements[i].nodes[face_id[1][0]],elements[i].nodes[face_id[1][1]],elements[i].nodes[face_id[1][2]],elements[i].nodes[face_id[1][3]]],elements[i].getDomainIndex(), -1))
        faces.append(Face(i,2,-1,0,[elements[i].nodes[face_id[2][0]],elements[i].nodes[face_id[2][1]],elements[i].nodes[face_id[2][2]],elements[i].nodes[face_id[2][3]]],elements[i].getDomainIndex(), -1))
        faces.append(Face(i,3,-1,0,[elements[i].nodes[face_id[3][0]],elements[i].nodes[face_id[3][1]],elements[i].nodes[face_id[3][2]],elements[i].nodes[face_id[3][3]]],elements[i].getDomainIndex(), -1))
        faces.append(Face(i,4,-1,0,[elements[i].nodes[face_id[4][0]],elements[i].nodes[face_id[4][1]],elements[i].nodes[face_id[4][2]],elements[i].nodes[face_id[4][3]]],elements[i].getDomainIndex(), -1))
        faces.append(Face(i,5,-1,0,[elements[i].nodes[face_id[5][0]],elements[i].nodes[face_id[5][1]],elements[i].nodes[face_id[5][2]],elements[i].nodes[face_id[5][3]]],elements[i].getDomainIndex(), -1))
    return faces

def getCornerNodeIDs(indexednodes,maxidx):
    
    imax = maxidx[0]+1
    jmax = maxidx[1]+1
    kmax = maxidx[2]+1
	
    i_index = 0
    j_index = 0
    k_index = 0
    node_000_indexedID = jmax*imax*k_index + imax*j_index + i_index

    i_index = imax-1
    j_index = 0
    k_index = 0
    node_x00_indexedID = jmax*imax*k_index + imax*j_index + i_index

    i_index = 0
    j_index = jmax-1
    k_index = 0
    node_0y0_indexedID = jmax*imax*k_index + imax*j_index + i_index

    i_index = imax-1
    j_index = jmax-1
    k_index = 0
    node_xy0_indexedID = jmax*imax*k_index + imax*j_index + i_index

    i_index = 0
    j_index = 0
    k_index = kmax-1
    node_00z_indexedID = jmax*imax*k_index + imax*j_index + i_index

    i_index = imax-1
    j_index = 0
    k_index = kmax-1
    node_x0z_indexedID = jmax*imax*k_index + imax*j_index + i_index

    i_index = 0
    j_index = jmax-1
    k_index = kmax-1
    node_0yz_indexedID = jmax*imax*k_index + imax*j_index + i_index

    i_index = imax-1
    j_index = jmax-1
    k_index = kmax-1
    node_xyz_indexedID = jmax*imax*k_index + imax*j_index + i_index

    node_000 = indexednodes[node_000_indexedID].getOldNodeID()
    node_x00 = indexednodes[node_x00_indexedID].getOldNodeID()
    node_0y0 = indexednodes[node_0y0_indexedID].getOldNodeID()
    node_xy0 = indexednodes[node_xy0_indexedID].getOldNodeID()
    node_00z = indexednodes[node_00z_indexedID].getOldNodeID()
    node_x0z = indexednodes[node_x0z_indexedID].getOldNodeID()
    node_0yz = indexednodes[node_0yz_indexedID].getOldNodeID()
    node_xyz = indexednodes[node_xyz_indexedID].getOldNodeID()

    domain_nodeID = (node_000, node_x00, node_0y0, node_xy0, node_00z, node_x0z, node_0yz, node_xyz)

    return domain_nodeID

def createIndexedNodes(nodes,style):

    maxidx = sizethegrid(nodes,style)

    imax = maxidx[0]+1
    jmax = maxidx[1]+1
    kmax = maxidx[2]+1

    print 'after sort imax =' + str(imax)
    print 'after sort jmax =' + str(jmax)
    print 'after sort kmax =' + str(kmax)
    
    indexednodes = []

    for node in range(len(nodes)):
        indicelist = nodes[node].indices
        if (len(indicelist)!=0):
            i_index=indicelist[0]
            j_index=indicelist[1]
            k_index=indicelist[2]
            indexednodeID = jmax*imax*k_index + imax*j_index + i_index
            indexednodes.append(IndexedNodes(indexednodeID,nodes[node].indices,nodes[node].position,nodes[node].id))

    sortID(indexednodes)

    return  indexednodes, maxidx

def renameIndexedNodes(nodes,style):

    maxidx = sizethegrid(nodes,style)

    imax = maxidx[0]+1
    jmax = maxidx[1]+1
    kmax = maxidx[2]+1

    print 'after sort imax =' + str(imax)
    print 'after sort jmax =' + str(jmax)
    print 'after sort kmax =' + str(kmax)
    
    indexednodes = []

    for node in range(len(nodes)):
        indicelist = nodes[node].indices
        if (len(indicelist)!=0):
            i_index=indicelist[0]
            j_index=indicelist[1]
            k_index=indicelist[2]
            indexednodeID = jmax*imax*k_index + imax*j_index + i_index
            indexednodes.append(IndexedNodes(indexednodeID,nodes[node].indices,nodes[node].position,nodes[node].oldnodeid))

    sortID(indexednodes)

    return  indexednodes, maxidx

def indexBlockSearch(domains, nnID,cellID,oldnodecellposition,nodes,nr,idv,idvo,style,celldelta,oldnodeID,domainID,previousDomainID,nodevisited,oddslice,oldnocells,sharedcells,firstblockpath,firstnodeflag):

    """
    Overview:    
    
    The local node zero is passed in from the outside.  The basic scheme is to set 
    the index of the current node, then craft vectors that specify how to find the 
    next node.  The number of neighboring cells surrounding each node in the block 
    is monitored.  Changes in the number of neighbors are used to sense when 
    transitions to/from center from/to face/edge/corner occur.  This information is 
    used to determine how to navigate through the block.  Each indexed node has a 
    flag set that it has been visited.  When all the nodes in the block are indexed, 
    then the sort terminates.
    
    Detail:        
        
    The process for the indexing can be broken down into indexing the node, then 
    finding the necessary vectors.  defineFaceandCellPosition determines for the
    current node (specified by nnID) where on the list of faces the node will be
    found.  This is fpv (face position vector), which is a 3-vector (for 3-D), each 
    position of which will contain a number between 0 and 5, corresponding to the 
    local x, y, and z directions.  nIDp (nodeID position) is a 4-vector (for 3D) 
    that specifies where in facenodelist nnID will be found, where each position has 
    a value between 0 and 3, (for 3D).  cellIDposition specifies SOMETHING.  
    
    The nfv (next face vector) is a simple opposite face restatement of fpv.  The nnv 
    (next node vector) is determined by recognizing that the next node will be on the 
    opposite face (specified by nfv) and will be on the nIDp for that face.  The nnv 
    is again a 3-vector, but each position contains the actual node ID corresponding 
    to the next node for each xyz direction.  
    
    nnID is moved to olnnodeID, and the idvUpdate is called.  idvUpdate is the main
    navigation routine within the indexblocksearch.  As mentioned in the overview, 
    the technique employed is that the number of cell neighbors surrounding the node
    is monitored.  If the number of cell neighbors has changed (embodied by the 
    parameter cellDelta) the center/face/edge/corner position in the block can be
    inferred.  The 3-vector idv (increment/decrement vector) specifies which direction
    will be next, and how the node will be indexed.  Once idv is established the nnID
    (new next Node) will be extracted from nnv.
    
    After nnID has been established, the number of neighboring (shared) cells and the 
    cellDelta is determined. The nodeVisited (i.e. indexed) flag is set.  Finally, the 
    previous (old) node cell position is established, along with the current cellID.
    The cellID will be used on the next pass in defineFaceandCellPosition.  If idvUpdate
    completed on the next pass (specified via the parameter done) then the for loop 
    terminates.
    """

    oldcellID = cellID
    firstnodeflag              = True

    for node in range(len(nodes)):

        nr                         = setIndices(nnID,cellID,oldnodecellposition,nodes,firstnodeflag,nr,idv,idvo,style,celldelta,oldnodeID,domainID)
        fpv, nIDp, cellIDposition  = defineFaceandCellPosition(domains, nodes, nnID, cellID, domainID, previousDomainID, firstblockpath, style)
        previousDomainID           = domainID
        nfv                        = buildNfv(fpv,style)
        nnv                        = buildNnv(nfv,nIDp,nodes,nnID,cellID,cellIDposition,domainID)
        oldnodeID                  = nnID
        idv, idvo,oddslice,done    = idvUpdate(idv,idvo,nnID,style,celldelta,nodevisited,oddslice,nnv,nodes)	   
        nnID                       = getNextNode(nnv,idv)

        firstnodeflag              = False

        oldnocells                 = sharedcells 
        celldelta, sharedcells     = cellDelta(oldnocells,nodes,nnID,domainID)
        nodevisited                = nodeVisited(nodes,nnID)
        if (cellID!=-1):
            oldcellID              = cellID

        oldnodecellposition,cellID = findCell(nnID,  oldcellID, nfv, nodes, idv,style, celldelta, nodevisited, domainID, firstblockpath)

        if done:
            break

        firstnodeflag              = False


def indexBlock(domains, nodes,domainID,previousDomainID,cornernodeID,cellID,style,nr,firstblockpath):
    
    """
    indexBlock sets up the calls to index the current block (specified by domainID).
    indexBlockSearch performs the actual indexing on the current block via the nodes
    object.  createIndexedNodes then extracts the just-indexed block of nodes from
    the object nodes, as a separate object (indexednodes) and returns said object 
    along with its size (maxidx).
    """

    idv = [1,0,0]
    idvo = idv

    firstnodeflag = False
#    firstendofrow = True
#    lengthofi  = 0
    oldnocells = 1
    oldnodeID  = cornernodeID
    oddslice= True

    oldnodecellposition = 0
    nnID = cornernodeID
    firstnodeflag = True
    celldelta = 0
    nodevisited = False
    sharedcells = 0
    previousDomainID = -1

    indexBlockSearch(domains, nnID,cellID,oldnodecellposition,nodes,nr,idv,idvo,style,celldelta,oldnodeID,domainID,previousDomainID, nodevisited,oddslice,oldnocells,sharedcells,firstblockpath,firstnodeflag)

    indexednodes, maxidx = createIndexedNodes(nodes,style)

    return  indexednodes, maxidx

def findSubsequentBlockCornersandFaces(domains, nodes, nextdomainIDlist, currentdomain, currentDomainID, face, style, firstblockpath):
    """ 
    Overview:    
    
    The purpose of findSubsequentBlockCornersandFaces is to find the corner nodes of
    the next block.  This next block is specified via domains[nextDomainID].  After
    determining the corners, the connectivity of the block with other blocks is 
    established.  The cornernodeID and cellID of node zero are established, before
    updating the previousDomainID.
    
    Detail:
        
    The indentification of this block is first determined as the nextDomainID.  The
    nodelist corresponding to the face of the block that is shared with the previous
    block (currentdomain).  This is the blocksharedfacenodelist.  Now in the nextDomain,
    the opposite face will be available.  To find this blockoppositefacenodelist, 
    buildOppositeFace is called.  With these two faces, whose node corners (n000, etc.)
    are exactly know, a geometry variable can be established.  This geometry variable
    lists out the corners in an exact specified order.  All this work is conducted in
    calculateDomainCorners.  
    
    After storing geometry in the domains object, the geometry is passed to 
    calculateFacesFromGeometry.  Here the face connectivity is established under the 
    domains[nextDomainID].blockfaceneighborlist.  The blockfacelist is also calculated
    in calculateFacesFromGeometry.  From here, retrieval of the cornernodeID is trivial.
    
    If there weren't any domains to actually inspect, i.e. len(nextdomainIDlist)==0, 
    then an error value for all the pertinent parameters is assigned.
    """
#    print 'findSubsequentBlockCornersandFaces:                      face = '+str(face)

    if (len(nextdomainIDlist)!=0):
        nextDomainID = nextdomainIDlist[0]
        blocksharedfacenodelist = currentdomain.blockfacelist[face]
        blockoppositefacenodelist = buildOppositeFace(domains,blocksharedfacenodelist,nextDomainID,nodes,style, firstblockpath)
        geometry = calculateDomainCorners(style, blocksharedfacenodelist, blockoppositefacenodelist, face,domains,nextDomainID,currentDomainID)
        domains[nextDomainID].geometry = geometry
        domains[nextDomainID].calculateFacesFromGeometry(style, geometry)#            cornernodeID = blocksharedfacenodelist[0]
            
        domainface_xneg = domains[nextDomainID].blockfacelist[1]
        cornernodeID = domainface_xneg[0]
        
        node = nodes[cornernodeID]
        for block in  range(len(node.domain)):
            blockID = node.domain[block]
            if (blockID==nextDomainID):
                cellID = node.neighbors[block]

        previousDomainID = currentDomainID
           
    else:
        print 'findSubsequentBlockCornersandFaces:  node domain not found on face '+str(face)+' of domain '+str(currentDomainID)
        previousDomainID = -1
        nextDomainID = -1
        cellID = -1
        cornernodeID = -1

    return  previousDomainID, nextDomainID, cellID, cornernodeID 


def recursiveCall(currentdomainID, nodes, domains,style, nr, filename):

    """
    This routine loops over the faces of a block (aka domain), checking to see if 
    there is another block connected to the face. If a block is attached to that 
    face, then the code determines if it has been indexed.  
    
    If the attached block has not been indexed, then the 
    findSubsequentBlockCornersandFaces establishes the corner nodes for that block, 
    returning the node and cell ZERO for that block.  It also returns the domanID's 
    as appropriate.
    
    If the cornernodeID is valid, this attached block is indexed using indexBlock.  
    This function returns the unrotated nodearray and maxidx.  maxidx is a 3 vector
    that contains the size of the nodearray.

    The unrotated nodearray is then subjected to the necessary rotations to match the
    indices to the previous block (domain).  The maxidx and rotatednodearray are 
    returned by the performRotations routine.

    The rotatednodearray is sent to the writeASCIIFile, where an individual file is 
    written for that block.  The next/previous domainID's are incremented, and the 
    recursiveCall is made.  
    
    When a block, during the recursiveCall, completes the scan of the faces conducted
    under the "for" loop, then routine returns to the previous block, which is also
    examining each of its faces.  That "for" loop also iterates toward completion.
    In this fashion, the block (domain) connectivity of the entire grid is explored.
    """

#    print 'recursiveCall: entering'
    currentdomain = domains[currentdomainID]
#    print 'recursiveCall: currentdomain = '+str(currentdomainID)
#    print 'recursiveCall: domains['+str(currentdomainID)+'].blockfacelist = ' +str(domains[currentdomainID].blockfacelist)
    for face in range(6):

        print 'recursiveCall: current face is '+str(face)

        nextdomainlistlength = len(currentdomain.blockfaceneighborlist[face])
        facepresentflag      = (nextdomainlistlength!=0)
        blockindexedflag     = False

        if (facepresentflag):
            nextdomainIDlist     = currentdomain.blockfaceneighborlist[face]
            nextdomain = domains[nextdomainIDlist[0]]
            blockindexedflag = nextdomain.blockindexed

        if (facepresentflag and not(blockindexedflag)):
            firstblockpath = True
            previousDomainID, nextDomainID, cellID, cornernodeID = findSubsequentBlockCornersandFaces(domains,nodes,nextdomainIDlist,currentdomain,currentdomainID, face,style, firstblockpath)

            if (cornernodeID != -1):

                nodearray, maxidx = indexBlock(domains, nodes,nextDomainID,previousDomainID, cornernodeID, cellID, style,nr, firstblockpath)

                rotatednodearray, maxidx = performRotations(nodearray,maxidx,nodes,domains,nextDomainID,style)

                print "recursiveCall: Saving Block..." +str(nextDomainID)        
                outfilename = filename + "_block_" +str(nextDomainID)+ ".inp"    
                writeASCIIFile(outfilename, rotatednodearray, nodes, style)
                print "recursiveCall: Written to " + outfilename
                currentdomainID = nextDomainID

                recursiveCall(currentdomainID,nodes,domains, style, nr, filename)
            else:
                print 'recursiveCall: face = '+str(face) +', does not have  a another block attached to it.'
        else:
            if (not(facepresentflag)):
                print 'recursiveCall: for domain '+str(currentdomainID)+ ', no domain detected on face = '+str(face)
            elif (blockindexedflag):
                print 'recursiveCall: for domain '+str(currentdomainID)+ ', domain '+str(nextdomain.id)+' detected on face = '+str(face)+' previously indexed'
    
    print 'recursiveCall: returning from domain  = '+str(currentdomainID)
    print 'recursiveCall: domainfaceneighborlist = '+str(domains[currentdomainID].blockfaceneighborlist)

def findNodeandFacePositionFromFacelist(facelist,nodeID):

    nodeID_nodeposvec = [0]*3
    nodeID_faceposvec = [0]*3

    facepos = 0

    for faceposition in range(len(facelist)):
        face = facelist[faceposition]
        nodeposition = returnNodePositionFromFaceNodeList(face,nodeID)
        if (nodeposition != -1):
            nodeID_nodeposvec[facepos] = nodeposition
            nodeID_faceposvec[facepos] = faceposition
            facepos = facepos + 1

    return nodeID_nodeposvec, nodeID_faceposvec
            
            
def find4thNodeinFace(face,nID):

    for nodeposition in range(len(face)):
        currentnodeID = face[nodeposition]
        test = (currentnodeID == nID[0])or(currentnodeID == nID[1])or(currentnodeID == nID[2])
        if (not(test)):
            nID_4th = currentnodeID
            break

    return nID_4th

def setIndices(nodeID,cellID,cellposition,nodes,firstnodeflag,nr,idv,idvo,style,deltacells,oldnodeID,domainID):
	
    """
    setIndices does the actual setting of the node indices.  If it is the first node, 
    then the value of the node indices is set based on an initialized len(idv)-vector
    If not the first node, then the idv value is applied to the node index of the 
    previous (old) node index to create the index for the current node (nodeID, also
    referred to external to setIndices as nnID).  The node indices are then added to 
    the node in question (i.e.nodes[nodeID]).  The number of nodes remaining is then
    decremented and returned as the value nr.  nr is no longer accurately intialized,
    is not used, and should be disregarded.
    """

    nodeindices =[0]*len(idv)

    if (firstnodeflag==True):
	    	print 'initialized nodeindices = ' + str(nodeindices)
    else:
		oldnodeindices = nodes[oldnodeID].indices
		for direction in range(len(nodeindices)):
		  nodeindices[direction] = oldnodeindices[direction]+idv[direction]  

    nodes[nodeID].addIndices(nodeindices)

    nnr = nr - 1

    return nnr 

def findFirstBlockCornerNode(domains, nodes,nr,style,domainID,firstblockpath):
	cornerfound = False
        
        position = getStartCornerFromKeyboard()
        refxpos = round(position[0],1)
        refypos = round(position[1],1)
        refzpos = round(position[2],1)

        print 'position from keyboard = ('+str(refxpos)+', '+str(refypos)+', '+str(refzpos)+')'

	for node in nodes:
#            singledomainneighborlist = node.findSharedCellsForSingleDomain(domainID)
#            singledomaingeometrylist = node.findSharedCellsGeometryForSingleDomain(domainID)
            totalneighborlist = node.getNeighbors()
            numofneighbors = len(totalneighborlist)

            """
            print 'findCornerNode:  totalneighborlist = ' +str(totalneighborlist)
            print 'findCornerNode:  numoftotalneighborlist = '+str(numofneighbors)
            print 'findCornerNode:  singledomainneighborlist = '+str(totalneighborlist)
            """

            if (numofneighbors==1):
                print 'findCornerNode:  node '+str(node.id)+' position = '+str(node.position)
                xpos = round(node.position[0],1)
                ypos = round(node.position[1],1)
                zpos = round(node.position[2],1)
                if ((xpos ==  refxpos) and (ypos == refypos) and (zpos == refzpos)):
                    print 'findCornerNode:  which is the cornernode we desire.'
                    cornerfound = True	    
                    cellID = totalneighborlist[0]
                    cornernodeID = node.id
                    geometrylist=node.geometry[0]
                    firstnodeflag=True
                    idv        = [1,0,0]
                    idvo       = [1,0,0]
                    deltacells = 0
                    oldnodeID = -1

                    domainID = node.domain[0]
                    nr       = setIndices(cornernodeID,cellID,geometrylist,nodes,firstnodeflag,nr,idv,idvo,style,deltacells,oldnodeID,domainID)
                    """
                    print 'findCornerNode:  cornernodeID   = '+str(cornernodeID)
                    print 'findCornerNode:  cellID         = '+str(cellID)
                    print 'findCornerNode:  domainID       = '+str(domainID)
                    print 'findCornerNode:  firstblockpath = '+str(firstblockpath)
                    print 'findCornerNode:  node['+str(cornernodeID)+'].domain = '+str(node.domain)
                    """
                    nfv, nnv      = setNfvandNnv(domains, nodes,cornernodeID,cellID,style,domainID,firstblockpath)
                    nextnodeID    = getNextNode(nnv,idv)
                    oldnodeID     = cornernodeID
                    firstnodeflag = False

                    break
 #               else:
 #                   print 'findCornerNode:  Cornernode found, but not what we want.  Continuing to the next node'
 #           else:
 #               print 'findCornerNode:  Cornernode not found.  Continuing to the next node'
	if (cornerfound==True):
	     print "findCornerNode:  desired corner nodeID is " + str(cornernodeID) + " found in cell " +str(cellID)+' in domain '+str(node.domain)
	else:
	     print "findCornerNode:  corner node was not found."
	return nextnodeID, cellID, nnv, idv, nfv, nr, cornernodeID, domainID 

def getStartCornerFromKeyboard():
    
    position = [0,0,0]

    for i in range(3):
        if (i == 0):
            posdir = 'x'
        elif (i == 1):
            posdir = 'y'
        elif (i == 2):
            posdir = 'z'
        nb = raw_input('enter the '+posdir+' position value for the first corner node:')

        try:
            position[i] = float(nb)
        except:
            print 'invalid number'
            
    return position

def establishBlindFPV(facelist,nodeID):
    
    """
    establishBlindFPV is nothing more than a search for nodeID
    in the facelist (recorded as fpv) and node position on that
    face (recorded as nIDp).  faces is the face nodelist.
    """

    fpv  = []
    nIDp = []
    faces= []

    for face in range(len(facelist)):
         currentface =facelist[face]
         for nodeposition in range(len(currentface)):
             nodefoundonface = (currentface[nodeposition]==nodeID)
             if nodefoundonface: 
                 fpv.append(face)
                 nIDp.append(nodeposition)
                 faces.append(currentface)

    return fpv, nIDp, faces

def   findFaces(previous_faces, current_faces, nodeID):

      numofpreviousfaces = len(previous_faces)
      numofcurrentfaces = len(current_faces)

      if ((numofpreviousfaces == 0)or(numofcurrentfaces==0)):
          print 'findFaces:  Error!  faces are missing!'

#      print 'findFaces: previousfaces = '+str(previous_faces)
#      print 'findFaces: currentfaces  = '+str(current_faces)


      for previous_face_position in range(numofpreviousfaces):
          previousfacenodelist = previous_faces[previous_face_position]
#          print 'findFaces: previousfacenodelist ='+str(previous_faces[previous_face_position])
          for current_face_position in range(numofcurrentfaces):
              currentfacenodelist = current_faces[current_face_position]
#              print 'findFaces: currentfacenodelist ='+str(current_faces[current_face_position])
              facematch = 0
              for previous_node_position in range(len(previousfacenodelist)):
                  previous_node = previousfacenodelist[previous_node_position]
                  for current_node_position in range(len(currentfacenodelist)):
                      current_node = currentfacenodelist[current_node_position]
#                      print 'findFaces:  previous_node, current_node = '+str(previous_node)+', '+str(current_node)
                      nodematch = (previous_node==current_node)
                      if (nodematch):
                          facematch = facematch+1
#                          print 'findFaces:  currentnode: facematch = '+str(current_node)+':  '+str(facematch)
                          if (facematch >= 3):
                               break
                  if (facematch >=3):
                      break
              if (facematch >=3):
                  break
#              print 'findFaces:  zeroing facematch'
              facematch = 0
          if (facematch >=3):
              break

      return previous_face_position, current_face_position, previous_node_position, current_node_position

def getOppositeFaceFromFacelist(faces,current_fpv,current_fp):
    
    shared_fp = current_fpv[current_fp]
    
    if (shared_fp==0):
        opposite_fp = 1
    elif (shared_fp==1):
        opposite_fp = 0
    elif (shared_fp==2):
        opposite_fp = 3
    elif (shared_fp==3):
        opposite_fp = 2
    elif (shared_fp==4):
        opposite_fp = 5
    elif (shared_fp==5):
        opposite_fp = 4

    """
    print 'getOppositeFaceFromFacelist:  faces = ' + str(faces)
    print 'getOppositeFaceFromFacelist:  shared_fp = '+str(shared_fp)
    print 'getOppositeFaceFromFacelist:  opposite_fp = '+str(opposite_fp)
    print 'getOppositeFaceFromFacelist:  shared_facenodelist = '+str(faces[shared_fp])
    print 'getOppositeFaceFromFacelist:  opposite_facenodelist = '+str(faces[opposite_fp])
    """
    oppositefacenodelist = faces[opposite_fp]

    return oppositefacenodelist

def rotateSharedFaceNodeList(oldcellsharedfacenodelist,newcellsharedfacenodelist):

      facelength = len(oldcellsharedfacenodelist)

      neg = [0, 0, 0, 0]

      for oldcellnodeposition in range(facelength):
          nodefound = False
          for newcellnodeposition in range(facelength):
              if (oldcellsharedfacenodelist[oldcellnodeposition]==newcellsharedfacenodelist[newcellnodeposition]):
                  neg[oldcellnodeposition]=oldcellsharedfacenodelist[oldcellnodeposition]
                  nodefound = True
          if (nodefound==False):
              print 'rotateSharedFaceNodeList:  shared faces do not match between cells'
              break
                  
      return neg

def cellHexahedronFaceIndices():

    face1nindices = (0,3,2,1)
    face2nindices = (0,4,7,3)
    face3nindices = (0,1,5,4)

    face1pindices = (4,7,6,5)
    face2pindices = (1,5,6,2)
    face3pindices = (3,2,6,7)

    return face1pindices,face1nindices,face2pindices,face2nindices,face3pindices,face3nindices

def rotateFaces(current_np, previous_np, current_fp, previous_fp, oldcellsharedfacenodelist, newcellsharedfacenodelist, newcelloppositefacenodelist):

#      facerotationvalue = current_np-previous_np
      facelength = len(oldcellsharedfacenodelist)
      pos = [0, 0, 0, 0]
      neg = [0, 0, 0, 0]
      correctedgeometrylist = [0, 0, 0, 0, 0, 0, 0, 0]


      xpos = (previous_fp==0)
      xneg = (previous_fp==1)
      ypos = (previous_fp==2)
      yneg = (previous_fp==3)
      zpos = (previous_fp==4)
      zneg = (previous_fp==5)


      """
      print 'rotateFaces:  oldcellsharedfacenodelist = '+ str(oldcellsharedfacenodelist)
      print 'rotateFaces:  newcelloppositefacenodelist = '+ str(newcelloppositefacenodelist)
      print 'rotateFaces:  newcellsharedfacenodelist = '+ str(newcellsharedfacenodelist)
      """

      for oldcellnodeposition in range(facelength):
          for newcellnodeposition in range(facelength):
              if (oldcellsharedfacenodelist[oldcellnodeposition]==newcellsharedfacenodelist[newcellnodeposition]):
                  neg[oldcellnodeposition]=newcellsharedfacenodelist[newcellnodeposition]
                  pos[oldcellnodeposition]=newcelloppositefacenodelist[newcellnodeposition]

      for newcellnodeposition1 in range(facelength):
          for newcellnodeposition2 in range(facelength):
              if (newcellnodeposition1!=newcellnodeposition2):
                  if (neg[newcellnodeposition1]==neg[newcellnodeposition2]):
#                      print 'rotateFaces:  no face match found during rotation, neg = '+str(neg)
#                      print 'rotateFaces:  no face match found during rotation, pos = '+str(pos)
                      neg = [0,0,0,0]
                      pos = [0,0,0,0]
              
#      print 'rotateFaces:  the finished ROTATED oppositefacenodelist = '+str(pos)
#      print 'rotateFaces:  the finished ROTATED sharedfacenodelist = '+str(neg)


      face1pindices,face1nindices,face2pindices,face2nindices,face3pindices,face3nindices = cellHexahedronFaceIndices()

      for i in range(4):
          if (xneg):
#              print 'rotateFaces:  xneg'# at i = '+str(i) + ', where'
#              print 'face1nindices['+str(i)+'] = '+str(face1nindices[i])+', and'
#              print 'face1pindices['+str(i)+'] = '+str(face1pindices[i])
              correctedgeometrylist[face1nindices[i]]=neg[i]
              correctedgeometrylist[face1pindices[i]]=pos[i]
 #             print 'rotateFaces:  correctedgeometrylist = '+str(correctedgeometrylist)
          elif (xpos):
 #             print 'rotateFaces:  xpos'
              correctedgeometrylist[face1nindices[i]]=neg[i]
              correctedgeometrylist[face1pindices[i]]=pos[i]
 #             print 'rotateFaces:  correctedgeometrylist = '+str(correctedgeometrylist)
          elif (yneg):
#              print 'rotateFaces:  yneg'
              correctedgeometrylist[face2nindices[i]]=neg[i]
              correctedgeometrylist[face2pindices[i]]=pos[i]
          elif (ypos):
#              print 'rotateFaces:  ypos'
              correctedgeometrylist[face2nindices[i]]=neg[i]
              correctedgeometrylist[face2pindices[i]]=pos[i]
          elif (zneg):
#              print 'rotateFaces:  zneg'
              correctedgeometrylist[face3nindices[i]]=neg[i]
              correctedgeometrylist[face3pindices[i]]=pos[i]
          elif (zpos):
#             print 'rotateFaces:  zpos'
              correctedgeometrylist[face3nindices[i]]=neg[i]
              correctedgeometrylist[face3pindices[i]]=pos[i]

#      print 'rotateFaces:  correctedgeometrylist = '+str(correctedgeometrylist)
      
      return correctedgeometrylist

def returnCellPositionFromNeighborList(neighborlist,cellID):

#    print 'returnCellPositionFromNeighborList: neighborlist = ' +str(neighborlist)
#    print 'returnCellPositionFromNeighborList: cellID       = ' +str(cellID)

    lengthofneighborlist = len(neighborlist)

    if (lengthofneighborlist == 0):
        print 'returnCellPositionfromNeighborList:  The numober of neighbors is ZERO.'
    else:
        for cellposition in range(lengthofneighborlist):
            neighbor = neighborlist[cellposition]
            if (neighbor==cellID):
                cellIDposition = cellposition
                break
                
    return cellIDposition

def returnNodePositionFromFaceNodeList(facenodelist,nodeID):

#    print 'returnCellPositionFromNeighborList: neighborlist = ' +str(neighborlist)
#    print 'returnCellPositionFromNeighborList: cellID       = ' +str(cellID)

    lengthoffacenodelist = len(facenodelist)
    nodeIDposition = -1

    if (lengthoffacenodelist == 0):
        print 'returnNodePositionfromFaceNodeList:  The numober of nodes on this face is ZERO.'
    else:
        for nodeposition in range(lengthoffacenodelist):
            currentnodeID = facenodelist[nodeposition]
            if (currentnodeID==nodeID):
                nodeIDposition = nodeposition
                break
                
    return nodeIDposition

def defineFaceandCellPosition(domains, nodes,nodeID,cellID,domainID,previousDomainID,firstblockflag,style):

    """
    This module establishes the fpv, nIDp, and cellIDpositon.  This is accomplished
    by first extracting the neighborlist and facelist corresponding to a single
    domain.  The total neighborlist (neighborlist) is extracted directly from the 
    nodes object.  The facelist is extracted from the single domain facelist once
    the single domain cell position has been determined.  With facelist, and the nodeID
    the fpv, nIDp, and faces are generated by establishBlindFPV.  faces is not used.  
    The desired values are returned.  
    """

    cellIDposition = -1
    singledomainneighborlist = nodes[nodeID].findSharedCellsForSingleDomain(domainID)
    singledomainfacelist     = nodes[nodeID].findSharedCellsFaceListForSingleDomain(domainID)

    neighborlist = nodes[nodeID].neighbors

    singledomaincellIDposition = returnCellPositionFromNeighborList(singledomainneighborlist,cellID)
    cellIDposition             = returnCellPositionFromNeighborList(neighborlist,cellID)

    facelist = singledomainfacelist[singledomaincellIDposition]

    fpv  = []
    nIDp = []

    if (firstblockflag):
        fpv, nIDp, faces = establishBlindFPV(facelist,nodeID)
	
    return fpv, nIDp, cellIDposition

def setNfvandNnv(domains, nodes,nodeID,cellID,style,domainID,firstblockpath):
#################################################################
#  the module sets the next face vector (nfv) and the next node 
#  vector (nnv).  It does this by identifying each face that has
#  the current node (nodeID) for a given cell (cellID).  Since the 
#  faces in for each cell are paired, and since the pairs are sequenced
#  in a list such that, for a hexahedron, with faces 1-6, the pairing 
#  goes as 1-2, 3-4, 5-6.  Similarly for the rectangle grid. The next
#  node for a given direction will be in the same position in the next
#  face vector as the original node was for the original face.
#################################################################

#    print 'setNfvandNnv: entering'

    previousDomainID = -1 

#    print 'setNfvandNnv: domainID = '+str(domainID)
#    print 'setNfvandNnv: entering defineFaceandCellPosition'
    fpv, nIDp, cellIDposition = defineFaceandCellPosition(domains, nodes,nodeID,cellID,domainID,previousDomainID,firstblockpath,style)
#    print 'setNfvandNnv: exiting defineFaceandCellPosition'
    
    #        print 'fpv = '+str(fpv)
    #        print 'nIDp = '+str(nIDp)
    #        print 'cellIDposition = '+str(cellIDposition)
    
    nfv = buildNfv(fpv,style)
    
    #        print 'from within setNfvandNnv, and after buildnfv call, nfv = '+str(nfv)
    
    nnv = buildNnv(nfv,nIDp,nodes,nodeID,cellID,cellIDposition,domainID)

#    print 'setNfvandNnv: exiting'

    return nfv, nnv


def buildNnv(nfv,nIDp,nodes,nodeID,cellID,cellIDposition,domainID):

	nnv = []

	for dir in range(len(nfv)):
		nfp = nfv[dir]
		nodeIDposition = nIDp[dir]
#                if (domainID==2):
 #                   print 'buildNnv:  nodeIDposition = '+str(nodeIDposition)
		nextnode = findNode(nfv,nfp,nodes,nodeID,cellID,cellIDposition,nodeIDposition,domainID)
		nnv.append(nextnode)

	return nnv

def buildNfv(fpv,style):

#    print 'from within buildNfv: fpv = '+str(fpv)

    if (style=='Rectangle'):
		if (fpv[0]==0)and(fpv[1]==2):
		  nfv=[1,3]
		elif (fpv[0]==0)and(fpv[1]==3):
		  nfv=[1,2]
		elif (fpv[0]==1)and(fpv[1]==2):
		  nfv=[0,1]
		elif (fpv[0]==1)and(fpv[1]==3):
		  nfv=[0,2]
    elif (style=='Hexahedron'):
		if ((fpv[0]==0)and(fpv[1]==2)and(fpv[2]==4)):
		  nfv=[1,3,5]
		elif ((fpv[0]==0)and(fpv[1]==2)and(fpv[2]==5)):
		  nfv=[1,3,4]
		elif ((fpv[0]==0)and(fpv[1]==3)and(fpv[2]==4)):
		  nfv=[1,2,5]
		elif ((fpv[0]==0)and(fpv[1]==3)and(fpv[2]==5)):
		  nfv=[1,2,4]
		elif ((fpv[0]==1)and(fpv[1]==2)and(fpv[2]==4)):
		  nfv=[0,3,5]
		elif ((fpv[0]==1)and(fpv[1]==2)and(fpv[2]==5)):
		  nfv=[0,3,4]
		elif ((fpv[0]==1)and(fpv[1]==3)and(fpv[2]==4)):
		  nfv=[0,2,5]
		elif ((fpv[0]==1)and(fpv[1]==3)and(fpv[2]==5)):
		  nfv=[0,2,4]

        # loop over all faces

#    print 'from withing buildNfv: nfv = '+str(nfv)

    return nfv

def findNode(nfv,faceposition,nodes,nodeID,cellID,cellIDposition,nodeIDposition,domainID):

#    sharedCellListFaces = nodes[nodeID].findSharedCellsFaceListForSingleDomain(domainID)

#    facelist = sharedCellListFaces[cellIDposition]# <<<<<-------THIS NEEDS TO BE CHECKED from here for multiblock

    facelist = nodes[nodeID].faces[cellIDposition]

#    print 'findNode:  facelist = '+str(facelist)

    facenodelist = facelist[faceposition]
    
    nextnode = facenodelist[nodeIDposition]

    """
    print 'findNode:  domainID \t= '+str(domainID)
    print 'findNode:  nodeID \t= '+str(nodeID)
    print 'findNode:  cellID \t= '+str(cellID)
    print 'findNode:  cellIDposition = '+str(cellIDposition)
    print 'findNode:  nodeIDposition = '+str(nodeIDposition)
    print 'findNode:  faceposition = '+str(faceposition)
    print 'findNode:  sharedCellListFaces = '+str(sharedCellListFaces)
    print 'findNode:  facelist = '+str(facelist)
    print 'findNode:  facenodelist = '+str(facenodelist)
    print 'findNode:  nextnode = ' +str(nextnode)
    """
    return nextnode

def endOfRow(oldnocells,nodes,nodeID):
	newnocells=len(nodes[nodeID].neighbors)
	delta=newnocells-oldnocells
	if delta < 0:
	   startnewrow=True
	else:
	   startnewrow=False
	return startnewrow, newnocells

def cellDelta(oldnocells,nodes,nodeID,domainID):
	newnocells= nodes[nodeID].findNumofSharedCells(domainID)
	delta=newnocells-oldnocells
	return delta, newnocells

def nodeVisited(nodes,nodeID):
	visited=(nodes[nodeID].indices!=[])
	return visited

def nodeUnVisited(nodes,nodeID):
	nodes[nodeID].indices=[]
	return 

def findDirection(idv):

	for dir in range(len(idv)):
		if (idv[dir]!=0):
			direction=dir

	return direction	


def findOldCell(nodes,nodeID,oldcellID,oldfaceposition,direction,domainID,firstblockpath): 

        singledomainneighbors = nodes[nodeID].findSharedCellsForSingleDomain(domainID)
        sharedcellsfacelist   = nodes[nodeID].findSharedCellsFaceListForSingleDomain(domainID)
#        facelist = nodes[nodeID].faces

        """
        print 'findOldCell:  nodeID                = ' +str(nodeID)
        print 'findOldCell:  singledomainneighbors = ' +str(singledomainneighbors)
        print 'findOldCell:  oldcellID             = '+str(oldcellID)
        print 'findOldCell:  sharedcellsfacelist   = '+str(sharedcellsfacelist)
        print 'findOldCell:  facelist              = '+str(facelist)
        """
	numberofneighbors = len(singledomainneighbors)
        
 #       print 'findOldCell:  numberofneighbors = '+str(numberofneighbors)

        cellposition = -1

	for cell in range(numberofneighbors):
            cellID = singledomainneighbors[cell]
            if (numberofneighbors==1):
                cellposition=0

#            print "findOldCell:  (cellID, oldcellID) = ("+str(cellID)+","+str(oldcellID)+")"

            cellmatch = (cellID==oldcellID)

            if (cellmatch):
                cellposition = cell
                oldfacelist  = sharedcellsfacelist[cell]
                oldface      = oldfacelist[oldfaceposition]
                break
#                print 'cellposition = '+str(cellposition)

        if (cellposition==-1):
            print 'findOldCell:  error! cell position unassigned'
            return
                        
#        print 'findOldCell: (cellposition, oldface, and oldfacelist) = ('+str(cellposition)+', '+str(oldface)+', '+str(oldfacelist)+')'
	return cellposition, oldface, oldfacelist

        
def findNewCell(nodes,nodeID,oldcellID,oldface,style):

	defaultcellID = -1
	
	defaultcellposition = -1

	for cell in range(len(nodes[nodeID].neighbors)):
	 	currentcellID = nodes[nodeID].neighbors[cell]
                if (len(nodes[nodeID].neighbors)==1):
                        cellposition  = cell
                        cellID = currentcellID
                        break
	 	celltest =(currentcellID!=oldcellID)
	 	if celltest:
		  facelist=nodes[nodeID].faces[cell]
		  cellID, cellposition = faceTest(style,oldface,facelist,nodes,nodeID,cell,defaultcellID,defaultcellposition)

		  if (cellID!=-1):
			break  
	
	return cellposition,cellID


def faceTest(style,oldface,facelist,nodes,nodeID,cell,defaultcellID,defaultcellposition,domainID):

	cellID = defaultcellID
	cellposition = defaultcellID

	for position in range(len(facelist)):
	   localface            =(facelist[position])
	   if (style=='Rectangle'):
		match1= ((localface[0]==oldface[0])and(localface[1]==oldface[1]))
		match2= ((localface[0]==oldface[1])and(localface[1]==oldface[0]))
		if (match1 or match2):
			cellID = nodes[nodeID].neighbors[cell]
			cellposition = cell
#			print 'Cell found! cellID = '+str(cellID)+' for nodeID = '+str(nodeID)+' with a cell position of '+str(cellposition)
			break

	   if (style=='Hexahedron'):
#                        print 'faceTest:  localface = '+str(localface)
#                        print 'faceTest:  oldface = '+str(oldface)
			match1= ((localface[0]==oldface[0])and(localface[1]==oldface[1])and(localface[2]==oldface[2])and(localface[3]==oldface[3]))
			match2= ((localface[1]==oldface[0])and(localface[2]==oldface[1])and(localface[3]==oldface[2])and(localface[0]==oldface[3]))
                        match3= ((localface[2]==oldface[0])and(localface[3]==oldface[1])and(localface[0]==oldface[2])and(localface[1]==oldface[3]))
                        match4= ((localface[3]==oldface[0])and(localface[0]==oldface[1])and(localface[1]==oldface[2])and(localface[2]==oldface[3]))
			if (match1 or match2 or match3 or match4):
                            singledomainneighbors =nodes[nodeID].findSharedCellsForSingleDomain(domainID)
                            cellID = singledomainneighbors[cell]
                            cellposition = cell
                            #			  print 'Cell found! cellID = '+str(cellID)+' for nodeID = '+str(nodeID)+' with a cell position of '+str(cellposition)
                            break

	return cellID, cellposition


def   findCell(nodeID,oldcellID,nfv,nodes,idv,style,deltacell,nodevisited,domainID,firstblockpath):
    
	"""
	So this node is important for the purpose of finding the cellID and cellposition.
      First we find the position of the old face in the facelist (i.e. oldfaceposition).
      This information is used to generate the cell position, oldface nodelist, and 
      oldfacelist.  If a transition to a face, or edge, or corner is reached, then the
      cellID will remain the same as the direction makes a u turn through the next three
      nodes.  If a u turn is not in process, i.e. if the oldcellID != cellID, then the 
      cellID and cellposition are extracted by faceTest.  An error value of -1 is 
      checked to see if faceTest succeeded.
	"""

	direction = findDirection(idv)
	oldfaceposition= nfv[direction]

	cellposition, oldface, oldfacelist = findOldCell(nodes,nodeID,oldcellID,oldfaceposition,direction,domainID,firstblockpath)

	cellID = -1
	
	defaultcellposition = -1

	singledomainneighbors = nodes[nodeID].findSharedCellsForSingleDomain(domainID)
	singledomainfaces     = nodes[nodeID].findSharedCellsFaceListForSingleDomain(domainID)

	cornercell    = (len(singledomainneighbors)==1)
	edgereached   = (deltacell<0)

	for cell in range(len(singledomainneighbors)):
	 	currentcellID = singledomainneighbors[cell]

	 	if (edgereached or cornercell): 
                    if (currentcellID == oldcellID):
                        cellposition  = cell
                        cellID = currentcellID
                 
	 	celltest =(currentcellID!=oldcellID)
	 	if celltest:
                    facelist=singledomainfaces[cell]

                    cellID, cellposition = faceTest(style,oldface,facelist,nodes,nodeID,cell,cellID,defaultcellposition,domainID)
                    if (cellID!=-1):
                        break  
        if (cellID==-1):
            print "FindCell was unsuccessful."

        return  cellposition, cellID

def findFaceList(nodes,nodeID,cellID,domainID):

    singledomainneighborlist = nodes[nodeID].findSharedCellsForSingleDomain(domainID)
    singledomainfacelist     = nodes[nodeID].findSharedCellsFaceListForSingleDomain(domainID)

    cellIDposition = -1

    for cellposition in range(len(singledomainneighborlist)):
        if (singledomainneighborlist[cellposition]==cellID):
            cellIDposition=cellposition
            break

    if (cellIDposition == -1):
        print "cell not found in findFaceList"

    facelist = singledomainfacelist[cellIDposition]

    return facelist

def twodIdvUpdate(idv,idvo,deltacells,oddslice):

    """
    This module detects if the end of a row has been reached, and sets the new 
    increment/decrement vector (nidv) appropriately.  The variable deltacells 
    indicates that the number of neighboring cells for the node under consideration 
    has dropped, indicating that an edge has been reached. Under such conditions it 
    calls for the idv to specify an increment not in the i-th direction, but in the 
    j-th direction.

    After that one step in the j-th direction, the idv is set to move in the negative 
    (or positive) i-th direction, based on history (idvo is the old idv from the 
    previous node.
    """

    nidv = [0, 0, 0]

    if (deltacells<0):
	   if (oddslice):
              idv1=1
	   else:
	      idv1=-1

    sliceshift = (idvo[2]!=1)

    edge2corner = (deltacells==-1)
    face2edge = (deltacells==-2)
    interior2face = (deltacells==-4)
    exit_transition = (edge2corner or face2edge or interior2face) 

    corner2edge   = (deltacells==1)
    edge2face     = (deltacells==2)
    face2interior = (deltacells==4)
    entrance_transition = (corner2edge or edge2face or face2interior)

    movinginx = (idv[0]!=0)
    movinginy = ((idv[0]==0)and(idv[1]!=0))

    previouslymovinginNegX=(idvo[0]==-1)
    previouslymovinginPosX=(idvo[0]==+1)

    nochangeinsharedcells = (deltacells==0)

    if (exit_transition):
            if (sliceshift):
                if(movinginx): # reached an end of row in the forward direction.  Increment/decrement row index.
                    nidv[0]=0
                    nidv[1]=idv1
                elif(movinginy):      # previously incremented/decremented row index, check direction of previous row motion (idvo[0]), reverse row direction.
                    if (previouslymovinginNegX): # if previous row direction negative, now positive.  Shut off row incr/decr.
                        nidv[0]=+1
                        nidv[1]=0
                    if (previouslymovinginPosX): #if previous row direction positive, now negative.  Shut off row incr/decr.
                        nidv[0]=-1
                        nidv[1]=0              
            else: # no sliceshift
                if(movinginx): # reached an end of row in the forward direction.  Increment/decrement row index.
                    nidv[0]=0
                    nidv[1]=idv1
                elif(movinginy): # previously incremented/decremented row index, check direction of previous row motion (idvo[0]), reverse row direction.
                    if (previouslymovinginNegX): #if previous row direction negative, now positive.  Shut off row incr/decr.
                        nidv[0]=+1
                        nidv[1]=0
                    if (previouslymovinginPosX): #if previous row direction positive, now negative.  Shut off row incr/decr.
                        nidv[0]=-1
                        nidv[1]=0              
    elif (entrance_transition): # if changing status from the other direction, then
	  if(movinginy):  # having just done a row change, record everything
	    if (previouslymovinginNegX):
	      nidv[0]=+1
	      nidv[1]=0
	    elif (previouslymovinginPosX):
	      nidv[0]=-1
	      nidv[1]=0
	  else: # if no row change, then just keep plugging
	      nidv=idv
    elif (nochangeinsharedcells): 
	  if(movinginy): # if have just done a row increment/decrement, then advance all counters as appropriate
	    if (previouslymovinginNegX):
	      nidv[0]=+1
	      nidv[1]=0
	    if (previouslymovinginPosX):
	      nidv[0]=-1
	      nidv[1]=0
	  else: # if no row change, then just keep plugging
	      nidv=idv
    else:
	  nidv=idv
	
    return nidv 
###################################################################

def threedIdvUpdate(idv,idvo,nn,deltacells,visited,oddslice,nnv,nodes):

    """
    This is the 3d update to the increment/decrement vector.
    If the next node under consideration has been visited, a 3d (k-th)
    index increment is specified.  If not, then if a k-th index 
    increment was specified, then it is reset to allow a 2-d (i-th)
    index increment.  Control is then handed off to the 2dIdvUpdate
    routine to manage the 2d progression.
    """

    nidv=[0,0,0]
    
    newoddslice=oddslice

    movinginz = (idv[2]==+1)

    done = False

    if (visited):
	if (oddslice==True):
	  newoddslice=False
	else:
	  newoddslice=True
        nidv[2]=+1
        nidv[1]=0
        nidv[0]=0
	nidvo = idvo
    elif (movinginz):
        nidv[2]=0
        nidv[1]=0
        if (idvo[0]==+1):
	    nidv[0]=-1
        if (idvo[0]==-1):
	    nidv[0]=+1
	nidvo = idv
    elif (not(movinginz)):
        nidv=twodIdvUpdate(idv,idvo,deltacells,oddslice)
	nidvo = idv
    else:
        print "error in 3dIDVupdate!"

    nnID = getNextNode(nnv,nidv)
    nodevisited = nodeVisited(nodes,nnID)
    if (nodevisited):
        print "********************"
        print "nodevisited in idvUpdate = "+str(nodevisited)
        print "********************"
	if (oddslice==True):
	  newoddslice=False
	else:
	  newoddslice=True
        nidv[2]=+1
        nidv[1]=0
        nidv[0]=0
	nidvo = idvo
        nnID = getNextNode(nnv,nidv)
        nodevisited = nodeVisited(nodes,nnID)
        if (nodevisited):
            print "********************"
            print "idvUpdate Finished = "+str(nodevisited)
            print "********************"
            done = True
            return nidv, newoddslice, nidvo, done             


#    print "nidv(0,1,2) =("+str(nidv[0])+","+str(nidv[1])+","+str(nidv[2])+"),\t nidvo(0,1,2) =("+str(nidvo[0])+","+str(nidvo[1])+","+str(nidvo[2])+")"
    return nidv, newoddslice, nidvo, done 

def idvUpdate(idv,idvo,nn,style,deltacells,visited,oddslice,nnv,nodes):
    
    """    
    idvUpdate updates the increment/decrement vector (idv), which is used
    to control how the indexing proceeds.  This outer layer determines
    if the block is a 2-D or 3-D problem.  The 3-D  idvupdate calls the
    2-D idvupdate.
    """    
        
    if (style=='Rectangle'):
		nidv = twodIdvUpdate(idv,idvo,deltacells,oddslice)
    elif style=='Hexahedron':
		nidv,newoddslice, nidvo, done = threedIdvUpdate(idv,idvo,nn,deltacells,visited,oddslice,nnv,nodes)
                

    nidvo = idv

    return nidv, nidvo, newoddslice, done

def getNextNode(nnv,idv):
	nnfound = False
	nn = -1
	for pos in range(len(idv)):
	  if (idv[pos]!=0):
#	    print 'nnv='+str(nnv)
	    nn = nnv[pos]
	    nnfound = True
	    break
        if (nnfound==False):
	  print 'next node not found in nnv in getNextNode'
       
  	return nn

def writeASCIIFile(outfilename,indexednodes,nodes, style):

    # First we check if the file exists
    if os.access(outfilename, os.F_OK):
        os.remove( outfilename)

    # Next we create a new file
    f = open(outfilename, 'a')

    maxidx = sizethegrid(nodes,style)

    imax = maxidx[0]+1
    jmax = maxidx[1]+1
    kmax = maxidx[2]+1

    indexednodeIDmax = kmax*jmax*imax + kmax*jmax + kmax

    print 'indexednodeIDmax = '+str(indexednodeIDmax)
    print 'len(indexednodes) = '+str(len(indexednodes))


    for i in range(imax):
        for j in range(jmax):
            for k in range(kmax):
#                print '(i,j,k) = ('+str(i)+', '+str(j)+', '+str(k)+')'
                indexednodeID = kmax*jmax*i + kmax*j + k
 #               print 'indexednodeID = ' +str(indexednodeID)
                xposition = indexednodes[indexednodeID].position[0]
                yposition = indexednodes[indexednodeID].position[1]
                zposition = indexednodes[indexednodeID].position[2]
                indexnodeIDstring = str(indexednodes[indexednodeID].id)
                print str(indexednodes[indexednodeID].id)+'\t '+str(indexednodes[indexednodeID].oldnodeid)+'\t'+str(indexednodes[indexednodeID].indices)+'\t'+str(round(xposition,3))+',\t '+str(round(yposition,3))+',\t '+str(round(zposition,3))
                i_index = indexednodes[indexednodeID].indices[0]
                j_index = indexednodes[indexednodeID].indices[1]
                k_index = indexednodes[indexednodeID].indices[2]
                indexstring = str(i_index)+','+str(j_index)+','+str(k_index)
                positionstring = str(xposition)+','+str(yposition)+','+str(zposition)
                nodestring = indexnodeIDstring + ','+indexstring+','+positionstring+'\n'
 #               print nodestring
                f.write(nodestring)
 
    for node in range(len(nodes)):
        nodes[node].indices=[]

    f.close()



##################### Merging Sort Algorithm ######################
# This is how we map the connectivity of our grid
def sortFaces(faces):

    # This is a simple merge sort
    # It is not especially fast since it allocates a ton of memory
    # and moves data everywhere all the time, but I'll fix that later
        
    def merge_sort_faces(faces):
        if len(faces) < 2:
            return faces
        else:
            halfpoint = len(faces) / 2
            vec1 = merge_sort_faces(faces[:halfpoint])
            vec2 = merge_sort_faces(faces[halfpoint:])
            return face_merge(vec1,vec2)
        
    def face_merge(vec1,vec2):

        result = []
        i,j = 0,0
        l1 = len(vec1)
        l2 = len(vec2)
        while(i < l1 or j < l2):
            if(i < l1 and j < l2):
                choice = vec1[i].compareFaces(vec2[j])
                if choice == -1:
                    result.append(vec1[i])
                    i+=1
                elif choice == 1:
                    result.append(vec2[j])
                    j+=1
                elif choice == 0:
                    vec1[i].n2 = vec2[j].n1
                    vec1[i].n2_f = vec2[j].n1_f
                    vec1[i].n2_d = vec2[j].n1_d
                    result.append(vec1[i])
                    i+=1
                    j+=1
                else:
                    return result
            elif i >= l1:
                return result + vec2[j:]
            elif j >= l2:
                return result + vec1[i:]
            else:
                return result
        return result

    return merge_sort_faces(faces)

###########################################################
## Now we look into how to run through the system

# This is the algorithm that assigns a Global index to all of the Elements
def elementIndexer(elements, style):

    # Recursively running through these elements is going to require a massive recursion process
    # I have separated them out into Line, Rectangle and Hexahedral methods in order to streamline them as much as possible


    # This recursive algorithm is what runs through the Elements and assignes values
    def recursive3DElementIndexer(elementIndex, globalIndex, elements):
        if elements[elementIndex].setIndices(globalIndex):
            # -x side:
            if elements[elementIndex].neighbors[0] >= 0:
                globalIndex[0] += 1
                recursive3DElementIndexer(elements[elementIndex].neighbors[0], globalIndex, elements)
                globalIndex[0] -= 1
                
            # +x side:
            if elements[elementIndex].neighbors[0] >= 0:
                globalIndex[0] -= 1
                recursive3DElementIndexer(elements[elementIndex].neighbors[1], globalIndex, elements)
                globalIndex[0] += 1
                
            # -y side:
            if elements[elementIndex].neighbors[0] >= 0:
                globalIndex[1] += 1
                recursive3DElementIndexer(elements[elementIndex].neighbors[2], globalIndex, elements)
                globalIndex[1] -= 1
                
            # +y side:
            if elements[elementIndex].neighbors[0] >= 0:
                globalIndex[1] -= 1
                recursive3DElementIndexer(elements[elementIndex].neighbors[3], globalIndex, elements)
                globalIndex[1] += 1
                
            # -z side:
            if elements[elementIndex].neighbors[0] >= 0:
                globalIndex[2] += 1
                recursive3DElementIndexer(elements[elementIndex].neighbors[4], globalIndex, elements)
                globalIndex[2] -= 1
                
            # +z side:
            if elements[elementIndex].neighbors[0] >= 0:
                globalIndex[2] -= 1
                recursive3DElementIndexer(elements[elementIndex].neighbors[5], globalIndex, elements)
                globalIndex[2] += 1
    def recursive2DElementIndexer(elementIndex, globalIndex, elements):
        if elements[elementIndex].setIndices(globalIndex):
            # -x side:
            if elements[elementIndex].neighbors[0] >= 0:
                globalIndex[0] += 1
                recursive2DElementIndexer(elements[elementIndex].neighbors[0], globalIndex, elements)
                globalIndex[0] -= 1
                
            # +x side:
            if elements[elementIndex].neighbors[0] >= 0:
                globalIndex[0] -= 1
                recursive2DElementIndexer(elements[elementIndex].neighbors[1], globalIndex, elements)
                globalIndex[0] += 1
                
            # -y side:
            if elements[elementIndex].neighbors[0] >= 0:
                globalIndex[1] += 1
                recursive2DElementIndexer(elements[elementIndex].neighbors[2], globalIndex, elements)
                globalIndex[1] -= 1
                
            # +y side:
            if elements[elementIndex].neighbors[0] >= 0:
                globalIndex[1] -= 1
                recursive2DElementIndexer(elements[elementIndex].neighbors[3], globalIndex, elements)
                globalIndex[1] += 1

    def recursive1DElementIndexer(elementIndex, globalIndex, elements):
        if elements[elementIndex].setIndices(globalIndex):
            # -x side:
            if elements[elementIndex].neighbors[0] >= 0:
                globalIndex[0] += 1
                recursive1DElementIndexer(elements[elementIndex].neighbors[0], globalIndex, elements)
                globalIndex[0] -= 1
                
            # +x side:
            if elements[elementIndex].neighbors[0] >= 0:
                globalIndex[0] -= 1
                recursive1DElementIndexer(elements[elementIndex].neighbors[1], globalIndex, elements)
                globalIndex[0] += 1
    
    # We start with a random Element
    startingElementIndex = random.randint(0, len(elements))
    
    if style == "Line":
        recursive1DElementIndexer(startingElementIndex,[0],elements)
    elif style == "Rectangle":
        recursive2DElementIndexer(startingElementIndex,[0,0],elements)
    elif style == "Hexahedron":
        recursive3DElementIndexer(startingElementIndex,[0,0,0],elements)

def sizethegrid(nodes,style):
    currentimax = -1
    currentjmax = -1
    currentkmax = -1

    for node in range(len(nodes)):
        index = nodes[node].indices
        dim = len(index)
        if (dim!=0):
            if style == "Line":
                imaxtest = (currentimax < index[0])
                if imaxtest:
                    currentimax=index[0]
                maxidx = [currentimax]
            elif style == "Rectangle":
                imaxtest = (currentimax < index[0])
                if imaxtest:
                    currentimax=index[0]
                jmaxtest = (currentjmax < index[1])
                if jmaxtest:
                    currentjmax=index[1]
                maxidx = [currentimax,currentjmax]
            elif style == "Hexahedron":
#                print style
                imaxtest = (currentimax < index[0])
                if imaxtest:
                    currentimax=index[0]
                jmaxtest = (currentjmax < index[1])
                if jmaxtest:
                    currentjmax=index[2]
                kmaxtest = (currentkmax < index[2])
                if kmaxtest:
                    currentkmax=index[2]
                maxidx = [currentimax,currentjmax,currentkmax]
            else:
                print spacing + "This is an unknown primitive for a structured grid"



    print "maxidx =" + str(maxidx)

    return maxidx



# First we take our Grid and map out the connectivity
def constructNodes(style, nodes, elements, domains, filename):

    spacing = "\t"
    
    print spacing + "Transferring element neighbor data to node list ..."
    
    time0 = time()
    
#    faces = []

    """
    """
    for domain in range(len(domains)):
	print 'domain identify =' +str(domains[domain].getID())
 

#    print 'constructNodes: a nuisance print statement'

    for element in elements:
	for nodeid in element.nodes:
	    nodes[nodeid].addNeighbor(element.id,element.nodes,style,element.domain,domains)
            

#    for domain in range(len(domains)):
#        for node in nodes:
#            exteriordomain = domains[domain].getID()
#            numofsharedcells = node.findNumofSharedCells(exteriordomain)

    time1 = time()

    print spacing + "Element neighbor data transfered to node list in " +str(time1-time0)+"seconds"

    nr = len(nodes)

    print "There are " + str(nr) + " nodes to begin."
    
    domainID = 0

    print "The current block is " +str(domainID)

    print 'constructNodes:  entering findCornerNode'
    firstblockpath = True
    nnID, cellID, nnv, idv, nfv, nr, cornernodeID, domainID = findFirstBlockCornerNode(domains, nodes,nr,style,domainID,firstblockpath)
    previousDomainID = -1
    print 'constructNodes:  exiting findCornerNode'

    nodearray, maxidx = indexBlock(domains, nodes,domainID,previousDomainID,cornernodeID,cellID,style,nr,firstblockpath)

    domains[domainID].generateBlockCornerNodes(nodearray, maxidx, style)
    print 'constructNodes:  domains[domainID].geometry with original corner nodes  = '+str(domains[domainID].geometry)

    adjustedgeometry = adjustCornerNodesForFirstBlock(domains,nodes,domainID,nodearray)
    domains[domainID].geometry = adjustedgeometry
    print 'constructNodes:  domains[domainID].geometry with adjusted corner nodes  = '+str(domains[domainID].geometry)
    domains[domainID].calculateFacesFromGeometry(style, adjustedgeometry)
    
    print 'constructNodes:  blockfacelist = '+str(domains[domainID].blockfacelist)

    rotatednodearray, maxidx = performRotations(nodearray,maxidx,nodes,domains,domainID,style)

    print 'constructNodes:  blockfaceneighborlist = '+str(domains[domainID].blockfaceneighborlist)

    print "Saving Block..." +str(domainID)        
    outfilename = filename + "_block_" +str(domainID)+ ".inp"    
    writeASCIIFile(outfilename,rotatednodearray, nodes, style)
        
    recursiveCall(domainID,nodes,domains, style, nr,filename)

    time1 = time()

    print spacing + "Nodes indexed in " + str(time1 - time0) + " seconds"
    
    return 



def adjustCornerNodesForFirstBlock(domains,nodes,domainID,nodearray):

    n0 =domains[domainID].n000
    n1 =domains[domainID].nx00
    n2 =domains[domainID].n0y0
    n3 =domains[domainID].n00z
    n4 =domains[domainID].nxy0
    n5 =domains[domainID].n0yz
    n6 =domains[domainID].nx0z
    n7 =domains[domainID].nxyz

    cornernodes = [n0,n1,n2,n3,n4,n5,n6,n7]

    maxnumberofnodes = len(cornernodes)

    npyna1D = zeros(maxnumberofnodes*3)

    ndist = npyna1D.reshape(maxnumberofnodes,3)

    for j in range(8):
        for i in range(3):
            ndist[j][i] = nodes[cornernodes[j]].position[i]

    ndtx =ndist[:,0]        
    ndty =ndist[:,1]        
    ndtz =ndist[:,2]        

    xminvalue=ndtx.min(0)
    yminvalue=ndty.min(0)
    zminvalue=ndtz.min(0)
    minxcorners = []
    minycorners = []
    minzcorners = []
        
    xmaxvalue=ndtx.max(0)
    ymaxvalue=ndty.max(0)
    zmaxvalue=ndtz.max(0)
    maxxcorners = []
    maxycorners = []
    maxzcorners = []
        
    for corner in range(len(ndtx)):
        if (ndtx[corner]==xminvalue):
            minxcorners.append(corner)
        if (ndty[corner]==yminvalue):
            minycorners.append(corner)
        if (ndtz[corner]==zminvalue):
            minzcorners.append(corner)
        if (ndtx[corner]==xmaxvalue):
            maxxcorners.append(corner)
        if (ndty[corner]==ymaxvalue):
            maxycorners.append(corner)
        if (ndtz[corner]==zmaxvalue):
            maxzcorners.append(corner)

    """
    print 'minxcorners = '+str(minxcorners)
    print 'minycorners = '+str(minycorners)
    print 'minzcorners = '+str(minzcorners)

    print 'maxxcorners = '+str(maxxcorners)
    print 'maxycorners = '+str(maxycorners)
    print 'maxzcorners = '+str(maxzcorners)
    """

    for xc_min in minxcorners:
        for yc_min in minycorners:
            for zc_min in minzcorners:
                for xc_max in maxxcorners:
                    for yc_max in maxycorners:
                        for zc_max in maxzcorners:
                            n000test=((xc_min==yc_min)and(yc_min==zc_min))
                            nx00test=((xc_max==yc_min)and(yc_min==zc_min))
                            n0y0test=((xc_min==yc_max)and(yc_max==zc_min))
                            n00ztest=((xc_min==yc_min)and(yc_min==zc_max))
                            nxy0test=((xc_max==yc_max)and(yc_max==zc_min))
                            n0yztest=((xc_min==yc_max)and(yc_max==zc_max))
                            nx0ztest=((xc_max==yc_min)and(yc_min==zc_max))
                            nxyztest=((xc_max==yc_max)and(yc_max==zc_max))
                            if n000test:
                                n000_id = cornernodes[xc_min]
                            if nx00test:
                                nx00_id = cornernodes[xc_max]
                            if n0y0test:
                                n0y0_id = cornernodes[yc_max]
                            if n00ztest:
                                n00z_id = cornernodes[zc_max]
                            if nxy0test:
                                nxy0_id = cornernodes[zc_min]
                            if n0yztest:
                                n0yz_id = cornernodes[xc_min]
                            if nx0ztest:
                                nx0z_id = cornernodes[xc_max]
                            if nxyztest:
                                nxyz_id = cornernodes[xc_max]

                                
    """
    print 'n000id = '+str(n000_id)
    print 'nx00id = '+str(nx00_id)
    print 'n0y0id = '+str(n0y0_id)
    print 'n00zid = '+str(n00z_id)
    print 'nxy0id = '+str(nxy0_id)
    print 'n0yzid = '+str(n0yz_id)
    print 'nx0zid = '+str(nx0z_id)
    print 'nxyzid = '+str(nxyz_id)

    print 'n000idpos = '+str(nodes[n000_id].position)
    print 'nx00idpos = '+str(nodes[nx00_id].position)
    print 'n0y0idpos = '+str(nodes[n0y0_id].position)
    print 'n00zidpos = '+str(nodes[n00z_id].position)
    print 'nxy0idpos = '+str(nodes[nxy0_id].position)
    print 'n0yzidpos = '+str(nodes[n0yz_id].position)
    print 'nx0zidpos = '+str(nodes[nx0z_id].position)
    print 'nxyzidpos = '+str(nodes[nxyz_id].position)
    """

    domains[domainID].addN000(n000_id)
    domains[domainID].addNx00(nx00_id)
    domains[domainID].addN0y0(n0y0_id)
    domains[domainID].addN00z(n00z_id)
    domains[domainID].addNxy0(nxy0_id)
    domains[domainID].addN0yz(n0yz_id)
    domains[domainID].addNx0z(nx0z_id)
    domains[domainID].addNxyz(nxyz_id)

    n000id =domains[domainID].n000
    nx00id =domains[domainID].nx00
    n0y0id =domains[domainID].n0y0
    n00zid =domains[domainID].n00z
    nxy0id =domains[domainID].nxy0
    n0yzid =domains[domainID].n0yz
    nx0zid =domains[domainID].nx0z
    nxyzid =domains[domainID].nxyz

    geometry = [n000id,n0y0id,n0yzid,n00zid,nx00id,nxy0id,nxyzid,nx0zid]

#    print 'cornernodes after processing = '+str(cornernodes)

    return geometry

def performRotations(nodearray,maxidx,nodes,domains,domainID,style):

    print 'performRotations:'
    domains[domainID].printDomainCornerNodePositions(nodes)

    numpynodearray           = translateNodeArrayToNumpy(nodearray,maxidx)
    rotatednumpynodearray    = rotateNumpyNodeArray(nodes,domains,numpynodearray,domainID)
    rotatednodearrayoldid    = translateNumpyNodeArrayBackToNodeArray(nodes,rotatednumpynodearray,nodearray,maxidx)
    rotatednodearray, maxidx = renameIndexedNodes(rotatednodearrayoldid,style)
    numpynodearray_newid     = translateNodeArrayToNumpy(rotatednodearray,maxidx)

    domains[domainID].blockindexed = True

    domains[domainID].findDomainNeighbor(nodes)
    domains[domainID].printDomainCornerNodeIDs(nodes)
    domains[domainID].printDomainCornerGeometryNodeIDs(nodes)

    print "performRotations: Storing indexed block in domain ..." +str(domainID)        
    domains[domainID].indexedBlock = numpynodearray_newid
    domains[domainID].maxidx = maxidx
    domains[domainID].getAndStoreIndexedBlockFaceSlices()


    return rotatednodearray, maxidx



def translateNumpyNodeArrayBackToNodeArray(nodes,rotatednumpynodearray,nodearray,maxidx):

    imax = maxidx[0]+1
    jmax = maxidx[1]+1
    kmax = maxidx[2]+1

    rotatednodearray = deepcopy(nodearray)

    for i in range(imax):
        for j in range(jmax):
            for k in range(kmax):
                oldnodeid = rotatednumpynodearray[i][j][k].oldnodeid

                rotatednodearray[oldnodeid].indices[0]     = i
                rotatednodearray[oldnodeid].indices[1]     = j
                rotatednodearray[oldnodeid].indices[2]     = k

                rotatednodearray[oldnodeid].position[0]    = rotatednumpynodearray[i][j][k].xpos              
                rotatednodearray[oldnodeid].position[1]    = rotatednumpynodearray[i][j][k].ypos             
                rotatednodearray[oldnodeid].position[2]    = rotatednumpynodearray[i][j][k].zpos                 
                 
  
    for rotatednode in rotatednodearray:
        i = rotatednode.indices[0]
        j = rotatednode.indices[1]
        k = rotatednode.indices[2]
        newindexednodeID = kmax*jmax*i + kmax*j + k
        rotatednode.id = newindexednodeID

    sortID(rotatednodearray)

#    print 'from within translate from numpy:'
#    printIndexedBlocks(rotatednodearray,nodearray)

    return rotatednodearray

def printIndexedBlocks(rotatednodearray,nodearray):

    for rotatednode in range(len(rotatednodearray)):
        print 'printIndexedBlocks:  rotatednodearray['+str(rotatednode)+'].oldnodeid = '+str(rotatednodearray[rotatednode].oldnodeid)
    for rotatednode in range(len(rotatednodearray)):
        print 'printIndexedBlocks:         nodearray['+str(rotatednode)+'].oldnodeid = '+str(nodearray[rotatednode].oldnodeid)

    return

def nodeMatchTest(nodes,domains,domainID):

    n000id = domains[domainID].n000
    n0y0id = domains[domainID].n0y0
    n0yzid = domains[domainID].n0yz
    n00zid = domains[domainID].n00z

#    print 'n000id = '+str(n000id)

    n000 = nodes[n000id].indices
    n0y0 = nodes[n0y0id].indices
    n0yz = nodes[n0yzid].indices
    n00z = nodes[n00zid].indices

    n000_test_i =    (n000[0]==0)
    n000_test_j =    (n000[1]==0)
    n000_test_k =    (n000[2]==0)

    n0y0_test_i =    (n0y0[0]==0)
    n0y0_test_j = not(n0y0[1]==0)
    n0y0_test_k =    (n0y0[2]==0)

    n0yz_test_i =    (n0yz[0]==0)
    n0yz_test_j = not(n0yz[1]==0)
    n0yz_test_k = not(n0yz[2]==0)

    n00z_test_i =    (n00z[0]==0)
    n00z_test_j =    (n00z[1]==0)
    n00z_test_k = not(n00z[2]==0)

    n000_test = n000_test_i*pow(2,2)+n000_test_j*pow(2,1)+n000_test_k*pow(2,0)

    n0y0_test = n0y0_test_i*pow(2,2)+n0y0_test_j*pow(2,1)+n0y0_test_k*pow(2,0)

    n0yz_test = n0yz_test_i*pow(2,2)+n0yz_test_j*pow(2,1)+n0yz_test_k*pow(2,0)

    n00z_test = n00z_test_i*pow(2,2)+n00z_test_j*pow(2,1)+n00z_test_k*pow(2,0)

#    print "nodeMatchTest:  n000_test(i,j,k) = ("+str(n000_test_i)+', '+str(n000_test_j)+', '+str(n000_test_k)+')'
#    print "nodeMatchTest:  n0y0_test(i,j,k) = ("+str(n0y0_test_i)+', '+str(n0y0_test_j)+', '+str(n0y0_test_k)+')'
#    print "nodeMatchTest:  n0yz_test(i,j,k) = ("+str(n0yz_test_i)+', '+str(n0yz_test_j)+', '+str(n0yz_test_k)+')'
#    print "nodeMatchTest:  n00z_test(i,j,k) = ("+str(n00z_test_i)+', '+str(n00z_test_j)+', '+str(n00z_test_k)+')'


    xneg_test=n000_test*pow(10,3)+n0y0_test*pow(10,2)+n0yz_test*pow(10,1)+n00z_test*pow(10,0)

    return xneg_test

def rotateNumpyNodeArray(nodes,domains,numpynodearray,domainID):

    rotatednumpynodearray = numpynodearray.copy()

    dfm = nodeMatchTest(nodes,domains,domainID)

    print 'dfm = '+str(dfm)

    if (dfm == 7777):
        print 'nothing to be done for domain '+str(domainID)
        axis1 = 1
        axis2 = 1
        rotatednumpynodearray=numpynodearray.swapaxes(axis1,axis2)
    elif (dfm == 7474):
        axis1 = 1
        axis2 = 2
        rotatednumpynodearray=numpynodearray.swapaxes(axis1,axis2)
    elif (dfm == 7117):
        axis1 = 0
        axis2 = 1
        rotatednumpynodearray=numpynodearray.swapaxes(axis1,axis2)
    elif (dfm == 7412):
        axis1 = 0
        axis2 = 1
        rotatednumpynodearray=numpynodearray.swapaxes(axis1,axis2)
        axis1 = 1
        axis2 = 2
        rotatednumpynodearray=rotatednumpynodearray.swapaxes(axis1,axis2)
    elif (dfm == 7722):
        axis1 = 0
        axis2 = 2
        rotatednumpynodearray=numpynodearray.swapaxes(axis1,axis2)
    elif (dfm == 7124):
        axis1 = 2
        axis2 = 0
        rotatednumpynodearray=numpynodearray.swapaxes(axis1,axis2)
        axis1 = 2
        axis2 = 1
        rotatednumpynodearray=rotatednumpynodearray.swapaxes(axis1,axis2)
    else:
        print 'rotateNumpyNodeArray:  Error!  Error code does NOT match'

    """
    print 'numpynodearray = '
    print numpynodearray
    print 'rotatednumpynodearray = '
    print rotatednumpynodearray
    """
    return rotatednumpynodearray

def translateNodeArrayToNumpy(nodearray,maxidx):

    imax = maxidx[0]+1
    jmax = maxidx[1]+1
    kmax = maxidx[2]+1

#    imaxm1 = imax - 1
#    jmaxm1 = jmax - 1
#    kmaxm1 = kmax - 1

    numofparams = 5
    numofdims = 3
    maxnumberofnodes = len(nodearray)

    npyna1Dfloat = zeros(maxnumberofnodes,float)
    npyna1Dint = zeros(maxnumberofnodes,int)

    npyna3D_nodeID      = npyna1Dint.reshape(imax,jmax,kmax)
    npyna3D_xpos        = npyna1Dfloat.reshape(imax,jmax,kmax)
    npyna3D_ypos        = npyna1Dfloat.reshape(imax,jmax,kmax)
    npyna3D_zpos        = npyna1Dfloat.reshape(imax,jmax,kmax)
    npyna3D_oldnodeid   = npyna1Dint.reshape(imax,jmax,kmax)
    
    npyna3D = np.core.records.fromarrays([npyna3D_nodeID,npyna3D_xpos,npyna3D_ypos,npyna3D_zpos,npyna3D_oldnodeid],names='nodeID,xpos,ypos,zpos,oldnodeid')

    for i in range(imax):
        for j in range(jmax):
            for k in range(kmax):
                indexednodeID = kmax*jmax*i + kmax*j + k
                i_index = nodearray[indexednodeID].indices[0]
                j_index = nodearray[indexednodeID].indices[1]
                k_index = nodearray[indexednodeID].indices[2]
                nodeID = nodearray[indexednodeID].id
                npyna3D[i_index][j_index][k_index][0]=nodeID
                npyna3D[i_index][j_index][k_index][1]=nodearray[indexednodeID].position[0]
                npyna3D[i_index][j_index][k_index][2]=nodearray[indexednodeID].position[1]
                npyna3D[i_index][j_index][k_index][3]=nodearray[indexednodeID].position[2]
                npyna3D[i_index][j_index][k_index][4]=nodearray[indexednodeID].oldnodeid

#    print 'translateNA2NMPY:  numpynodearray3D = \n'
#    print npyna3D

    """
    for k in range(kmax):
        for j in range(jmax):
            for i in range(imax):
                nodeID = npyna3D[i,j,k,0]
                print 'npyna3D('+str(i)+','+str(j)+','+str(k)+') =['+str(nodeID)+']'
    """

    return npyna3D

def main(argv):

    ngc = 1
    filename = ""
    formatting = ""   
    # We need to parse in the readfile
    if len(argv) < 3 or len(argv) > 3:
        print """"I need an input terms. Please set:
    python gridconverterhdf5.py inputfilename.whatever [ngc-x,ngc+x,ngc-y...]

Where ngc = number of ghost cells
        """
        return
    else:
        period = argv[1].index(".")
        filename = argv[1][:period]
        formatting = argv[1][period:]
        ngc = list(eval(argv[2]))

    time0 = time()

    print "Loading Input File..."
    
    nodes, elements, domains, style = readFile(filename + formatting, ngc)

    time1 = time()
    
    print "Loaded in " + str(time1 - time0) + " seconds"

    if(style == ""):
        print "Error Found. Exiting Script..."
        return

    print "Constructing Node Indices..."

    constructNodes(style, nodes, elements, domains, filename)

    del elements

    time2 = time()

    print "Nodes Constructed in " + str(time2-time1) + " seconds"

    time3 = time()

    print "File Saved in " + str(time3 - time2) + " seconds"
    
    print "Grid Rebuild Complete in " + str(time3-time0) + " seconds."

        
if __name__ == "__main__":
	main(sys.argv)
