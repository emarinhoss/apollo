#!/usr/bin/env python
r"""Creates XMF files for use in Visit
"""

header = r"""<?xml version="1.0" ?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>
<Xdmf Version="2.0">
 <Domain>
   <Grid Name="mesh1" GridType="Uniform">
"""

footer = r"""
   </Grid>
 </Domain>
</Xdmf>
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
        %s:/%s/Y
       </DataItem>
       <DataItem Dimensions="%d" NumberType="Float" Precision="4" Format="HDF">
        %s:/%s/X
       </DataItem>
     </Geometry>
"""

geometry3d = r"""
     <Geometry GeometryType="X_Y_Z">
       <DataItem Dimensions="%d" NumberType="Float" Precision="4" Format="HDF">
        %s:/%s/X
       </DataItem>
       <DataItem Dimensions="%d" NumberType="Float" Precision="4" Format="HDF">
        %s:/%s/Y
       </DataItem>
       <DataItem Dimensions="%d" NumberType="Float" Precision="4" Format="HDF">
        %s:/%s/Z
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
    parser.add_option('-w', '--writer', action='store',
                      dest='writer',
                      help='Name of SubSolver which writes data')
    parser.add_option('-d', '--dump-numbers', action='store',
                      dest='dumps',
                      help='Dumps to process. For example: -d 2,3 will process dumps 2 and 3'
                      ' To specify a range use -d "range(0,10)". This will process dumps'
                      ' 0-9 inclusive.')

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

    # check if write block exists
    if options.writer not in top.childrenMap:
        print("Writer %s does not exist in file") % options.writer
        exit(1)
    
    # get writer block
    writer = top.childrenMap[options.writer]
    # now determine the names of the variables written by this block
    outNms = writer.values['writeNames']
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
        # write out header data
        xmfFile.writelines(header)

        empty = []
        myOutNms = eval(outNms[0])
        if type(myOutNms) == type(""):
            myOutNms = (myOutNms,)

        if ndim == 2:
            # write toplogy data
            xmfFile.writelines(topology2d % (spatialOrder*cells[0], spatialOrder*cells[1]))
            # write grid information
            h5 = '%s_%d.h5' % (inpFile[:-4], i)
            xmfFile.writelines(geometry2d % (spatialOrder*cells[1], h5, top.name,
                                             spatialOrder*cells[0], h5, top.name))

            # now write out each data array
            for nm in myOutNms:
                xmfFile.writelines(attribute2d % (nm, spatialOrder*cells[0], spatialOrder*cells[1],
                                                  h5, top.name, nm))
        elif ndim == 3:
            # write toplogy data
            xmfFile.writelines(topology3d % (spatialOrder*cells[0],
                                             spatialOrder*cells[1],
                                             spatialOrder*cells[2]))
            # write grid information
            h5 = '%s_%d.h5' % (inpFile[:-4], i)
            xmfFile.writelines(geometry3d % (spatialOrder*cells[2], h5, top.name,
                                             spatialOrder*cells[1], h5, top.name,
                                             spatialOrder*cells[0], h5, top.name))

            # now write out each data array
            for nm in myOutNms:
                xmfFile.writelines(attribute3d % (nm,
                                                  spatialOrder*cells[0],
                                                  spatialOrder*cells[1],
                                                  spatialOrder*cells[2],
                                                  h5, top.name, nm))

        # write out footer data
        xmfFile.writelines(footer)

        xmfFile.close()
