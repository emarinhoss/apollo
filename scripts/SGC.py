import os
import os.path
import sys
#import h5py
import numpy
import random
import pdb
from time import clock, time

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

###########################################        

    def addNeighbor(self,cell, nodes, style):
	self.neighbors.append(cell)
	self.geometry.append(list(nodes))

	geometrylist = list(nodes)

	if (style == 'Hexahedron'):
	  face1pindices = (0,3,2,1)
	  face1nindices = (4,7,6,5)
	  face2pindices = (3,2,6,7)
	  face2nindices = (0,1,5,4)
	  face3pindices = (1,5,6,2)
	  face3nindices = (0,4,7,3)


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

	self.faces.append(faces)

###########################################        

    def addIndices(self,indices):
	self.indices=list(indices)

class IndexedNodes(object):

    def __init__(self, ident, index, position):
        self.id      = ident            # decimal number corresponding to indices
        self.indices = list(index)      # Logical grid location of the node
        self.position = list(position)  # Physical position of Node

    def getID(self):
        return self.id

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
            domainIndex = 0   # This will be the ID of each domain through reconstruction
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
                        print("This grid converter does not work with the Element primitive in the input file")
                        return [],[],[],0

                    # Now we append this new Domain to our list of Domains
                    domains.append(Domain(domainIndex, domain_name, style))
                    domains[-1].setGhostCells(ghostCells)
                    domainIndex += 1
                    
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
    print(spacing + "Nodes Found: " + str(len(nodes)))
    print(spacing + "Domains Found: " + str(len(domains)))
    print(spacing + "Elements Found: " + str(len(elements)))

    # The last thing to do before we return this data is to
    # sort the Nodes by their ID and map the Elements to
    # the Node Indices instead of the Node IDs

    print(spacing + "Sorting Nodes...")
    sortID(nodes)

    print(spacing + "Nodes Sorted")
    print(spacing + "Reindexing Element Nodes...")
    for i in range(len(elements)):
        nodeset = elements[i].nodes
        for j in range(len(nodeset)):
            nodeIndex = getIndexFromID(nodes,nodeset[j])
            elements[i].nodes[j] = nodeIndex
        #print("Rebuilding Element ") + str(i) + ": " + str(elements[i].nodes)
    print(spacing + "Element Nodes Reindexed")
    
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

def setIndices(nodeID,cellID,cellposition,nodes,firstnodeflag,nr,idv,idvo,style,deltacells,oldnodeID):
	
#	print('entering setIndices, where idv =')+str(idv)

    	nodeindices =[0]*len(idv)
	if (firstnodeflag==True):
	    	print('initialized nodeindices = ') + str(nodeindices)
	else:
		oldnodeindices = nodes[oldnodeID].indices
		facelist = findFaceList(nodes,nodeID,cellID)
		for direction in range(len(nodeindices)):
		  nodeindices[direction] = oldnodeindices[direction]+idv[direction]  

	nodes[nodeID].addIndices(nodeindices)

	print('nodeID = ') + str(nodeID) + ', where nodeindices after assignment = '+str(nodes[nodeID].indices)

	nnr = nr - 1

	return nnr 

def findCornerNode(nodes,nr,style):
	cornerfound = False
	for node in nodes:
	    numofneighbors=len(node.neighbors)
	    if (numofneighbors<3):
		    print(numofneighbors)
	    if (numofneighbors==1):
		cornerfound = True	    
		cellID = node.neighbors[0]
		cornernodeID = node.id
		geometrylist=node.geometry[0]
		firstnodeflag=True
		idv  = [1,0,0]
		idvo = [1,0,0]
		deltacells = 0
		oldnodeID = -1

	        nr = setIndices(cornernodeID,cellID,geometrylist,nodes,firstnodeflag,nr,idv,idvo,style,deltacells,oldnodeID)

		nfv, nnv = setNfvandNnv(nodes,cornernodeID,cellID,style)
     		nextnodeID = getNextNode(nnv,idv)
		oldnodeID = cornernodeID


		firstnodeflag=False

		break

	if (cornerfound==True):
	     print("corner nodeID is ") + str(cornernodeID) + " found in cell " +str(cellID)
	else:
	     print("corner node was not found.")
	return nextnodeID, cellID, nnv, idv, nfv, nr, cornernodeID 

def defineFaceandCellPosition(nodes,nodeID,cellID):

	cellIDposition = -1

	for cellposition in range(len(nodes[nodeID].neighbors)):
		if (nodes[nodeID].neighbors[cellposition]==cellID):
			cellIDposition = cellposition

	facelist =nodes[nodeID].faces[cellIDposition]

	fpv  = []
	nIDp = []

	for face in range(len(facelist)):
	   currentface =facelist[face]
	   for nodeposition in range(len(currentface)):
		nodefoundonface = (currentface[nodeposition]==nodeID)
		if nodefoundonface: 
		   fpv.append(face)
		   nIDp.append(nodeposition)

	return fpv, nIDp, cellIDposition

def setNfvandNnv(nodes,nodeID,cellID,style):
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

	fpv, nIDp, cellIDposition = defineFaceandCellPosition(nodes,nodeID,cellID)
	
	nfv = buildNfv(fpv,style)

        nnv = buildNnv(nfv,nIDp,nodes,nodeID,cellID,cellIDposition)

	return nfv, nnv

def buildNnv(nfv,nIDp,nodes,nodeID,cellID,cellIDposition):

	nnv = []

	for dir in range(len(nfv)):
		nfp = nfv[dir]
		nodeIDposition = nIDp[dir]
		nextnode = findNode(nfv,nfp,nodes,nodeID,cellID,cellIDposition,nodeIDposition)
		nnv.append(nextnode)

	return nnv

def buildNfv(fpv,style):

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

	return nfv

def findNode(nfv,faceposition,nodes,nodeID,cellID,cellIDposition,nodeIDposition):

	facelist = nodes[nodeID].faces[cellIDposition]

	facenodelist = facelist[faceposition]

	nextnode = facenodelist[nodeIDposition]

	return nextnode

def endOfRow(oldnocells,nodes,nodeID):
	newnocells=len(nodes[nodeID].neighbors)
	delta=newnocells-oldnocells
	if delta < 0:
	   startnewrow=True
	else:
	   startnewrow=False
	return startnewrow, newnocells

def cellDelta(oldnocells,nodes,nodeID):
	newnocells=len(nodes[nodeID].neighbors)
	delta=newnocells-oldnocells
	return delta, newnocells

def nodeVisited(nodes,nodeID):
	visited=(nodes[nodeID].indices!=[])
	return visited

def findDirection(idv):

	for dir in range(len(idv)):
		if (idv[dir]!=0):
			direction=dir

	return direction	

def findOldCell(nodes,nodeID,oldcellID,oldfaceposition,direction): 

	numberofneighbors = len(nodes[nodeID].neighbors)
	for cell in range(numberofneighbors):
	 	cellID = nodes[nodeID].neighbors[cell]
		if (numberofneighbors==1):
		   cellposition=0
		if (cellID==oldcellID):
			cellposition      = cell
			oldfacelist       = nodes[nodeID].faces[cell]
			oldface 	  = oldfacelist[oldfaceposition]

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


def faceTest(style,oldface,facelist,nodes,nodeID,cell,defaultcellID,defaultcellposition):

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
#			print('Cell found! cellID = ')+str(cellID)+' for nodeID = '+str(nodeID)+' with a cell position of '+str(cellposition)
			break

	   if (style=='Hexahedron'):
			match1= ((localface[0]==oldface[0])and(localface[1]==oldface[1])and(localface[2]==oldface[2])and(localface[3]==oldface[3]))
			match2= ((localface[1]==oldface[0])and(localface[2]==oldface[1])and(localface[3]==oldface[2])and(localface[0]==oldface[3]))
                        match3= ((localface[2]==oldface[0])and(localface[3]==oldface[1])and(localface[0]==oldface[2])and(localface[1]==oldface[3]))
                        match4= ((localface[3]==oldface[0])and(localface[0]==oldface[1])and(localface[1]==oldface[2])and(localface[2]==oldface[3]))
			if (match1 or match2 or match3 or match4):
			  cellID = nodes[nodeID].neighbors[cell]
			  cellposition = cell
#			  print('Cell found! cellID = ')+str(cellID)+' for nodeID = '+str(nodeID)+' with a cell position of '+str(cellposition)
			  break

	return cellID, cellposition

 
                
def findCell(nodeID,oldcellID,nfv,nodes,idv,style,deltacell,nodevisited):

	direction = findDirection(idv)

	oldfaceposition= nfv[direction]
	cellposition, oldface, oldfacelist = findOldCell(nodes,nodeID,oldcellID,oldfaceposition,direction)

	cellID = -1
	
	defaultcellposition = -1

	for cell in range(len(nodes[nodeID].neighbors)):
	 	currentcellID = nodes[nodeID].neighbors[cell]
                cornercell =len(nodes[nodeID].neighbors)==1
                edgereached = (deltacell<0)

                if (edgereached or cornercell):
                    if (currentcellID == oldcellID):
                        cellposition  = cell
                        cellID = currentcellID
                        break
	 	celltest =(currentcellID!=oldcellID)
	 	if celltest:
			facelist=nodes[nodeID].faces[cell]

			cellID, cellposition = faceTest(style,oldface,facelist,nodes,nodeID,cell,cellID,defaultcellposition)

			if (cellID!=-1):
#			  print('after faceTest, the cellID = ')+str(cellID)
			  break  
        if (cellID==-1):
		print("FindCell was unsuccessful.")

        return  cellposition, cellID

def findFaceList(nodes,nodeID,cellID):
	for cellposition in range(len(nodes[nodeID].neighbors)):
	   if (nodes[nodeID].neighbors[cellposition]==cellID):
           	cellIDposition=cellposition

	facelist = nodes[nodeID].faces[cellIDposition]

	return facelist

def twodIdvUpdate(idv,idvo,deltacells,oddslice):
###################################################################
# This module detects if the end of a row has been reached,
# and sets the new increment/decrement vector (nidv) appropriately.
# The variable deltacells indicates that the number of neighboring
# cells for the node under consideration has dropped, indicating that
# an edge has been reached. Under such conditions it calls for the idv
# to specify an increment not in the i-th direction, but in the j-th 
# direction.
#
# After that one step in the j-th direction, the idv is set to move in the 
# negative (or positive) i-th direction, based on history (idvo is the old
# idv from the previous node.
###################################################################
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
        endofrowreached = (movinginx and exit_transition and sliceshift)

        previouslymovinginNegX=(idvo[0]==-1)
        previouslymovinginPosX=(idvo[0]==+1)

        nochangeinsharedcells = (deltacells==0)

#        print("exit_transition, entrance_transition, sliceshift, movinginx, movinginy = (") + str(exit_transition) + ", " + str(entrance_transition)  + ", " + str(sliceshift) + ", "+str(movinginx) +","+str(movinginy)+")"

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
	elif(entrance_transition): # if changing status from the other direction, then
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
###################################################################
# This is the 3d update to the increment/decrement vector.
# If the next node under consideration has been visited, a 3d (k-th)
# index increment is specified.  If not, then if a k-th index 
# increment was specified, then it is reset to allow a 2-d (i-th)
# index increment.  Control is then handed off to the 2dIdvUpdate
# routine to manage the 2d progression.
###################################################################
    nidv=[0,0,0]
    
    newoddslice=oddslice

    movinginz = (idv[2]==+1)

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
        print("error in 3dIDVupdate!")

    nnID = getNextNode(nnv,nidv)
    nodevisited = nodeVisited(nodes,nnID)
    if (nodevisited):
        print("********************")
        print("nodevisited in idvUpdate = ")+str(nodevisited)
        print("********************")
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
            print("********************")
            print("idvUpdate Finished = ")+str(nodevisited)
            print("********************")
            return nidv, newoddslice, nidvo             


#    print("nidv(0,1,2) =(")+str(nidv[0])+","+str(nidv[1])+","+str(nidv[2])+"),\t nidvo(0,1,2) =("+str(nidvo[0])+","+str(nidvo[1])+","+str(nidvo[2])+")"
    return nidv, newoddslice, nidvo 

def idvUpdate(idv,idvo,nn,style,deltacells,visited,oddslice,nnv,nodes):
	if (style=='Rectangle'):
		nidv = twodIdvUpdate(idv,idvo,deltacells,oddslice)
	elif style=='Hexahedron':
		nidv,newoddslice, nidvo = threedIdvUpdate(idv,idvo,nn,deltacells,visited,oddslice,nnv,nodes)
                

	nidvo = idv

	return nidv, nidvo, newoddslice

def getNextNode(nnv,idv):
	nnfound = False
	nn = -1
	for pos in range(len(idv)):
	  if (idv[pos]!=0):
#	    print('nnv=')+str(nnv)
	    nn = nnv[pos]
	    nnfound = True
	    break
        if (nnfound==False):
	  print('next node not found in nnv in getNextNode')
       
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

    for i in range(imax):
        for j in range(jmax):
            for k in range(kmax):
                indexednodeID = kmax*jmax*i + kmax*j + k
                print(str(indexednodes[indexednodeID].id)+'\t'+str(indexednodes[indexednodeID].indices)+'\t'+str(indexednodes[indexednodeID].position))
                indexnodeIDstring = str(indexednodes[indexednodeID].id)
                i_index = indexednodes[indexednodeID].indices[0]
                j_index = indexednodes[indexednodeID].indices[1]
                k_index = indexednodes[indexednodeID].indices[2]
                indexstring = str(i_index)+','+str(j_index)+','+str(k_index)
                xposition = indexednodes[indexednodeID].position[0]
                yposition = indexednodes[indexednodeID].position[1]
                zposition = indexednodes[indexednodeID].position[2]
                positionstring = str(xposition)+','+str(yposition)+','+str(zposition)
                nodestring = indexnodeIDstring + ','+indexstring+','+positionstring+'\n'
 #               print(nodestring)
                f.write(nodestring)
 

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
            print(spacing + "This is an unknown primitive for a structured grid")
        
    return maxidx



# First we take our Grid and map out the connectivity
def constructNodes(style, nodes, elements, domains):

    spacing = "\t"
    
    print(spacing + "Transferring element neighbor data to node list ...")
    
    time0 = time()
    
    faces = []

    for domain in range(len(domains)):
	print('domain identify =') +str(domains[domain].getID())
	print(str(domains[domain].startIndices))
    for element in elements:
	for nodeid in element.nodes:
	    nodes[nodeid].addNeighbor(element.id,element.nodes,style)
    

    time1 = time()

    print(spacing + "Element neighbor data transfered to node list in " +str(time1-time0)+"seconds")
    
    nr = len(nodes)

    print(spacing + "There are " + str(nr) + " nodes to begin.")

    nodeID = -1
    nnID, cellID, nnv, idv, nfv,nr, cornernodeID = findCornerNode(nodes,nr,style)
    idvo = idv

    firstnodeflag = False
    firstendofrow = True
    lengthofi = 0
    oldnocells=1
    oldnodeID = cornernodeID
    oddslice= True

    for node in range(len(nodes)-1):
	celldelta, sharedcells = cellDelta(oldnocells,nodes,nnID)

#        print("celldelta = ") +str(celldelta)

	nodevisited                = nodeVisited(nodes,nnID)
        celllist                   = nodes[nnID].neighbors
        oldcellID                  = cellID
        oldnodecellposition,cellID = findCell(nnID,oldcellID,nfv,nodes,idv,style,celldelta,nodevisited)
        nr                         = setIndices(nnID,cellID,oldnodecellposition,nodes,firstnodeflag,nr,idv,idvo,style,celldelta,oldnodeID)

	fpv, nIDp, cellIDposition  = defineFaceandCellPosition(nodes,nnID,cellID)
	nfv                        = buildNfv(fpv,style)
        nnv                        = buildNnv(nfv,nIDp,nodes,nnID,cellID,cellIDposition)
	oldnodeID                  = nnID
	idv, idvo,oddslice         = idvUpdate(idv,idvo,nnID,style,celldelta,nodevisited,oddslice,nnv,nodes)	    
     	nnID                       = getNextNode(nnv,idv)

        veryoldnocells             = oldnocells
        oldcelldelta               = celldelta
        oldnocells                 = sharedcells 

    print("nodes indexed")

    maxidx = sizethegrid(nodes,style)

    imax = maxidx[0]+1
    jmax = maxidx[1]+1
    kmax = maxidx[2]+1

    print('after sort imax =') + str(imax)
    print('after sort jmax =') + str(jmax)
    print('after sort kmax =') + str(kmax)
    
    indexednodes = []

    for node in range(len(nodes)):
        indicelist = nodes[node].indices
        i_index=indicelist[0]
        j_index=indicelist[1]
        k_index=indicelist[2]
        indexednodeID = jmax*imax*k_index + imax*j_index + i_index
        indexednodes.append(IndexedNodes(indexednodeID,nodes[node].indices,nodes[node].position))
#        print str(indexednodeID)+'\t'+str(indexednodes[node].indices)+'\t'+str(indexednodes[node].position)

    sortID(indexednodes)

    nodearray = indexednodes 
    """
    for indexednode in range(len(indexednodes)):
        print(str(indexednodes[indexednode].id)+'\t'+str(indexednodes[indexednode].indices)+'\t'+str(indexednodes[indexednode].position))
    """
    """
    for nodeindex in range(len(nodearray)):
        print(str(nodearray[nodeindex].id)+'\t'+str(nodearray[nodeindex].indices)+'\t'+str(nodearray[nodeindex].position))
    """
    time1 = time()

    print(spacing + "Nodes indexed in " + str(time1 - time0) + " seconds")
    
    return nodearray


def main(argv):

    ngc = 1
    filename = ""
    formatting = ""   
    # We need to parse in the readfile
    if len(argv) < 3 or len(argv) > 3:
        print("")""I need an input terms. Please set:
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

    print("Loading Input File...")
    
    nodes, elements, domains, style = readFile(filename + formatting, ngc)

    time1 = time()
    
    print("Loaded in ") + str(time1 - time0) + " seconds"

    if(style == ""):
        print("Error Found. Exiting Script...")
        return

    print("Constructing Node Indices...")

    nodearray = constructNodes(style, nodes, elements, domains)

    del elements

    time2 = time()

    print("Nodes Constructed in ") + str(time2-time1) + " seconds"

    print("Saving File...")
    
    outfilename = filename + "_converted.inp"
    
    writeASCIIFile(outfilename,nodearray, nodes, style)

    time3 = time()

    print("File Saved in ") + str(time3 - time2) + " seconds"
    
    print("Grid Rebuild Complete in ") + str(time3-time0) + " seconds."
    print("Written to ") + outfilename
        
if __name__ == "__main__":
	main(sys.argv)
