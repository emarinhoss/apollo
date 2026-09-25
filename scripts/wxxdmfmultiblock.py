#!/usr/bin/env python
r"""Creates XMF files for use in Visit
"""

header = r"""<?xml version="1.0" ?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>
<Xdmf Version="2.0">
 <Domain>
"""

footer = r"""
 </Domain>
</Xdmf>
"""

gridheader = r"""
   <Grid Name="mesh%s" GridType="Uniform">
"""

gridfooter = r"""
   </Grid>
"""

topology2d = r"""
   <Topology TopologyType="2DRectMesh" NumberOfElements="%d %d 1"/>
"""

topology3d = r"""
   <Topology TopologyType="3DRectMesh" NumberOfElements="%d %d %d 1"/>
"""

geometry2d = r"""
     <Geometry GeometryType="X_Y_Z">
       <DataItem Dimensions="%d" NumberType="Float" Precision="4" Format="HDF">
        %s:/%s/Y%s
       </DataItem>
       <DataItem Dimensions="%d" NumberType="Float" Precision="4" Format="HDF">
        %s:/%s/X%s
       </DataItem>
     </Geometry>
"""

geometry3d = r"""
     <Geometry GeometryType="X_Y_Z">
       <DataItem Dimensions="%d" NumberType="Float" Precision="4" Format="HDF">
        %s:/%s/X%s
       </DataItem>
       <DataItem Dimensions="%d" NumberType="Float" Precision="4" Format="HDF">
        %s:/%s/Y%s
       </DataItem>
       <DataItem Dimensions="%d" NumberType="Float" Precision="4" Format="HDF">
        %s:/%s/Z%s
       </DataItem>       
     </Geometry>
"""

attribute2d = r"""
     <Attribute Name="%s" AttributeType="Scalar" Center="Node">
       <DataItem Dimensions="%d %d 1" NumberType="Float" Precision="4" Format="HDF">
         %s:/%s/%s
       </DataItem>
     </Attribute>
"""

attribute3d = r"""
     <Attribute Name="%s" AttributeType="Scalar" Center="Node">
       <DataItem Dimensions="%d %d %d 1" NumberType="Float" Precision="4" Format="HDF">
         %s:/%s/%s
       </DataItem>
     </Attribute>
"""

import wxinputparser

from optparse import OptionParser

if __name__ == '__main__':
    # read command line options
    parser = OptionParser()
    parser.add_option('-i', '--input', action='store',
                      dest='inputFile',
                      help='Input file to preprocess')
    parser.add_option('-d', '--dump-numbers', action='store',
                      dest='dumps',
                      help='Dumps to process. For example: -d 2,3 will process dumps 2 and 3'
                      ' To specify a range use -d "range(0,10)". This will process dumps'
                      ' 0-9 inclusive.')
    parser.add_option('-s', '--datawriterStep', action='store',
                      dest='datawriterstep',
                      help='Name of SubSolverStep which writes data')
    parser.add_option('-g', '--gridwriterStep', action='store',
                      dest='gridwriterstep',
                      help='Name of SubSolverStep which writes grid data')

    (options, args) = parser.parse_args()

    # determine name of input file
    if options.inputFile:
        inpFile = options.inputFile
    else:
        if len(args) == 0:
            print("Must provide name of input file.")
            exit(1)
        else:
            inpFile = args[0]

    # open input file for parsing
    if inpFile[-4:] == '.pin':
        inpFile = inpFile[:-4] + '.inp'
    elif inpFile[-4:] == '.inp':
        pass
    else:
        inpFile = inpFile + '.inp'

    inp = wxinputparser.WxInpParse(inpFile)

    # top level WarpX block
    top = inp.dataSet.children[0]

    # get writer block
    datawriterstep = top.childrenMap[options.datawriterstep]
    # now determine the names of the variables written by this block
    outNmsWS = datawriterstep.values['SubSolvers']
    numWriters = len(outNmsWS)
#    print(numWriters)
#    print outNmsWS[0]


    for i in range(numWriters):
        # get writer block
        writer = top.childrenMap[outNmsWS[i]]
        # now determine the names of the variables written by this block
        outNms = writer.values['writeNames']
#	print outNms[0]
        grid = top.childrenMap[writer.values['OnGrid']]
#	print(grid.values)
    # get writer block
    gridwriterstep = top.childrenMap[options.gridwriterstep]
    # now determine the names of the variables written by this block
    outNmsGS = gridwriterstep.values['SubSolvers']
    numGrids = len(outNmsGS)
#    print(numGrids)
#    print(outNmsGS)
    for i in range(numGrids):
        # get writer block
        gridwriter = top.childrenMap[outNmsGS[i]]
        # now determine the names of the variables written by this block
        outGNms = gridwriter.values['writeNames']
#	print("outGNms = ")+str(outGNms)


    # check what kind it is
    kind = writer.values["Kind"]
    spatialOrder = 1
    meqn = 1
    if kind == 'exprWriterDG':
        meqn = writer.values["meqn"]
        spatialOrder = writer.values["spatialOrder"]

    grid = top.childrenMap[writer.values['OnGrid']]

    # now fetch the grid size
    cells = grid.values['Cells']
    # determine dimensions
    ndim = len(cells)

    # determine number of output files
    nout = top.values['Out']

    # determine which ranges to work on
    if options.dumps:
        dumps = eval(options.dumps)
        if type(dumps) == type(0):
            dumps = (dumps,)
    else:
        dumps = range(nout+1)

    # now loop, writing out each XMF file
    for i in dumps:
        xmfFile = open('%s_%d.xmf' % (inpFile[:-4], i), 'w')
        h5 = '%s_%d.h5' % (inpFile[:-4], i)
        # write out header data
        xmfFile.writelines(header)

        empty = []

        if ndim == 2:
            for nm in range(numGrids):
                writer = top.childrenMap[outNmsWS[nm]]
                # now determine the names of the variables written by this block
                outNms = writer.values['writeNames']
                grid = top.childrenMap[writer.values['OnGrid']]
                # now fetch the grid size
                cells = grid.values['Cells']
                # write grid data
                xmfFile.writelines(gridheader % (nm+1))
                # write toplogy data
                xmfFile.writelines(topology2d % (spatialOrder*cells[0], 
                                                 spatialOrder*cells[1]))
                # write geometry data
                xmfFile.writelines(geometry2d % (spatialOrder*cells[1], h5, top.name,nm+1,
                                                 spatialOrder*cells[0], h5, top.name,nm+1))
                # extract the writeNames from the excessive quotation marks
                myOutNms = eval(outNms[0])
                if type(myOutNms) == type(""):
                        myOutNms = (myOutNms,)
                # now write out each data array
                for j in range(len(myOutNms)):
                        xmfFile.writelines(attribute2d % (myOutNms[0], 
                                                          spatialOrder*cells[0], 
                                                          spatialOrder*cells[1],
                                                          h5, top.name, 
                                                          myOutNms[0]))
                xmfFile.writelines(gridfooter)

        elif ndim == 3:
            for nm in range(numGrids):
                writer = top.childrenMap[outNmsWS[nm]]
                # now determine the names of the variables written by this block
                outNms = writer.values['writeNames']
                grid = top.childrenMap[writer.values['OnGrid']]
                # now fetch the grid size
                cells = grid.values['Cells']
                # write grid data
                xmfFile.writelines(gridheader % (nm+1))
                # write toplogy data
                xmfFile.writelines(topology3d % (spatialOrder*cells[0], 
                                                 spatialOrder*cells[1], 
                                                 spatialOrder*cells[2]))
                # write geometry data
                xmfFile.writelines(geometry3d % (spatialOrder*cells[2], h5, top.name,nm+1,
                                                 spatialOrder*cells[1], h5, top.name,nm+1,
                                                 spatialOrder*cells[0], h5, top.name,nm+1))
                # extract the writeNames from the excessive quotation marks
                myOutNms = eval(outNms[0])
                if type(myOutNms) == type(""):
                        myOutNms = (myOutNms,)
                # now write out each data array
                for j in range(len(myOutNms)):
                        xmfFile.writelines(attribute3d % (myOutNms[0], 
                                                          spatialOrder*cells[0],
                                                          spatialOrder*cells[1],
                                                          spatialOrder*cells[2],
                                                          h5, top.name, 
                                                          myOutNms[0]))
                xmfFile.writelines(gridfooter)

        # write out footer data
        xmfFile.writelines(footer)

        xmfFile.close()
