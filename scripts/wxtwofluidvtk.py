#!/usr/bin/env python

import wxdata
import numpy

from optparse import OptionParser

def writeHeaderAndGrid(fd, q):
    fd.write('# vtk DataFile Version 2.0\n')
    fd.write('WarpX Two-Fluid data\n')    
    fd.write('ASCII\n')

    grd = q.grid # get hold of grid object
    nx, ny, nz, me = q.shape
    
    # write out grid: this is non-uniform rectangular grid
    fd.write('DATASET RECTILINEAR_GRID\n')
    fd.write('DIMENSIONS %d %d %d\n' % (nx, ny, nz))
    # write X coordinates
    fd.write('X_COORDINATES %d float\n' % nx)
    grd.X.tofile(fd, '\n')

    # write Y coordinates
    fd.write('\nY_COORDINATES %d float\n' % ny)
    grd.Y.tofile(fd, '\n')

    # write Z coordinates
    fd.write('\nZ_COORDINATES %d float\n' % nz)
    grd.Z.tofile(fd, '\n')

def writeScalar(fd, data, name):
    nx, ny, nz = data.shape
    ndata = nx*ny*nz # total number of elements to write
    
    # write out scalar data
    fd.write("\nPOINT_DATA %d\n" % ndata) # header
    fd.write("\nSCALARS %s float 1\n" % name) # scalar data
    fd.write("LOOKUP_TABLE default\n")
    data.transpose().tofile(fd, "\n") # actual data
    fd.write("\n")

def writeVector(fd, vx, vy, vz, name):
    nx, ny, nz = vx.shape
    ndata = nx*ny*nz # total number of elements to write

    # write out vector data
    fd.write("\nPOINT_DATA %d\n" % ndata) # header
    fd.write("VECTORS %s float\n" % name) # vector data

    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                fd.write("%g %g %g\n" % (vx[i,j,k], vy[i,j,k], vz[i,j,k]))
    fd.write("\n")

def convertTwoFluid2VTK(inpFile, frame, arrName, meqn):
    r"""convertTwoFluid2VTK(inpFile, frame, arrName, meqn)

    Converts H5 data from a file into a VTK file for visit plotting.
    """

    # open data file and read in DG array
    d = wxdata.WxData(inpFile, frame)
    q = d.readDG(arrName, meqn)

    inpFile = inpFile + "_" + str(frame)

    # write electron density
    print("Writing electron density....")
    fn = inpFile + "_elc_density.vtk"
    fd = open(fn, "w")
    writeHeaderAndGrid(fd, q)
    writeScalar(fd, q[:,:,:,0], "elc_density")
    fd.close()

    # write electron velocity
    print("Writing electron velocity....")
    fn = inpFile + "_elc_velocity.vtk"
    fd = open(fn, "w")
    writeHeaderAndGrid(fd, q)
    vx = q[:,:,:,1]/q[:,:,:,0]
    vy = q[:,:,:,2]/q[:,:,:,0]
    vz = q[:,:,:,3]/q[:,:,:,0]
    writeVector(fd, vx, vy, vz, "elc_velocity")

    # write electron pressure
    print("Writing electron pressure....")
    fn = inpFile + "_elc_pressure.vtk"
    fd = open(fn, "w")
    writeHeaderAndGrid(fd, q)
    # compute pressure first
    pr = (5.0/3.0-1)*(q[:,:,:,4]
                      - 0.5*(q[:,:,:,1]*q[:,:,:,1]
                             + q[:,:,:,2]*q[:,:,:,2]
                             + q[:,:,:,3]*q[:,:,:,3])/q[:,:,:,0])
    writeScalar(fd, pr, "elc_pressure")
    fd.close()    

    # write ion density
    print("Writing ion density....")
    fn = inpFile + "_ion_density.vtk"
    fd = open(fn, "w")
    writeHeaderAndGrid(fd, q)
    writeScalar(fd, q[:,:,:,5], "ion_density")
    fd.close()

    # write ion velocity
    print("Writing ion velocity....")
    fn = inpFile + "_ion_velocity.vtk"
    fd = open(fn, "w")
    writeHeaderAndGrid(fd, q)
    vx = q[:,:,:,6]/q[:,:,:,5]
    vy = q[:,:,:,7]/q[:,:,:,5]
    vz = q[:,:,:,8]/q[:,:,:,5]
    writeVector(fd, vx, vy, vz, "ion_velocity")

    # write ion pressure
    print("Writing ion pressure....")
    fn = inpFile + "_ion_pressure.vtk"
    fd = open(fn, "w")
    writeHeaderAndGrid(fd, q)
    # compute pressure first
    pr = (5.0/3.0-1)*(q[:,:,:,9]
                      - 0.5*(q[:,:,:,6]*q[:,:,:,6]
                             + q[:,:,:,7]*q[:,:,:,7]
                             + q[:,:,:,8]*q[:,:,:,8])/q[:,:,:,5])
    writeScalar(fd, pr, "ion_pressure")
    fd.close()

    # write electric field
    print("Writing electric field....")
    fn = inpFile + "_elec.vtk"
    fd = open(fn, "w")
    writeHeaderAndGrid(fd, q)    
    ex = q[:,:,:,10]
    ey = q[:,:,:,11]
    ez = q[:,:,:,12]
    writeVector(fd, ex, ey, ez, "electric_field")
    fd.close()

    # write electric field
    print("Writing magnetic field....")
    fn = inpFile + "_magn.vtk"
    fd = open(fn, "w")
    writeHeaderAndGrid(fd, q)    
    bx = q[:,:,:,13]
    by = q[:,:,:,14]
    bz = q[:,:,:,15]
    writeVector(fd, bx, by, bz, "magnetic_field")
    fd.close()

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option('-i', '--input',
                      dest='inpFile',
                      help = 'Input file (no extesion) to convert')
    parser.add_option('-f', '--frame',
                      dest='frame',
                      help = 'Frame number')
    parser.add_option('-a', '--array',
                      dest='arrName',
                      help='Name of array')
    parser.add_option('-m', '--num-equations',
                      dest = 'meqn',
                      help = 'Number of equations')

    # parse command line arguments
    (options, args) = parser.parse_args()

    inpFile = options.inpFile
    frame = int(options.frame)
    arrName = options.arrName
    meqn = int(options.meqn)

    # perform conversion
    convertTwoFluid2VTK(inpFile, frame, arrName, meqn)
