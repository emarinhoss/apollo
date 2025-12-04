r"""Provides an interface to read and plot WarpX hyperbolic solver
data. The numpy package is required. The code works with 1D/2D/3D
simulations.

To read data from a run with input file 'testrun.inp' use::

 >>> td = wxdata.WxData('testrun.inp')

This will read meta-data from the input file and the simulation log
file and initialize an object which can be used to read the actual data
into memory. The returned object ``td`` has some useful
information. To determine the dimension of the simulation::

 >>> td.dims
 2

To determine the grid used::

 >>> td.cells
 [128, 64]

To determine the coordinates of the lower-left corner::

 >>> td.lower
 [-12.800000000000001, -6.4000000000000004]


To determine the coordinates of the upper-right corner::
 >>> td.upper
 [12.800000000000001, 6.4000000000000004]

To get the x-coordinates of cell  centers, for example::

 >>> xc = td.axis[0] # xc is a Numeric array
 >>> xc.shape
 (128,)

similarly for Y and Z coordinates.

Once the ``WxData`` object is created one can read the data from a
particular frame as follows::

 >>> q5 = td.read(5) # read data from frame 5. q0 is a Numeric array
 >>> q5.shape
 (128, 64, 5)

For 1D simulations the returned array is of shape ``(nx, meqn)`` and
for 2D simulation is of shape ``(nx, ny, meqn)``, where ``meqn`` is
the number of equations in the system. The read data can be now
manipulated or plotted as desired."""

import numpy
import re

class WxData:
    r"""WxData(inpFile : string) -> WxData

    Provides an interface to read data from a WarpX hyperbolic solver
    simulation with input file ``inpFile``. The data is returned as
    numpy array."""

    def __init__(self, inp, extra=1):

        # this is a HACK!! Should be reading input file as XML and
        # parsing data out of it
        
        # open input file and read contents
        data = open(inp).read()
        # search for Run_name
        pr = re.compile(r'Run_name\s*=\s*(?P<rn>[a-zA-Z0-9_\-]*)')
        self.base = pr.search(data).group('rn')
        # search for block name
        prb = re.compile(r'Blocks\s*=\s*\[\s*(?P<ar>.*)\]')
        self.block = prb.search(data).group('ar')
        # search for cells
        prc = re.compile(r'Cells\s*=\s*(?P<ar>\[.*\])')
        self.cells = eval(prc.search(data).group('ar'))
        # search for lower
        prl = re.compile(r'Lower\s*=\s*(?P<ar>\[.*\])')
        self.lower = eval(prl.search(data).group('ar'))
        # search for upper
        pru = re.compile(r'Upper\s*=\s*(?P<ar>\[.*\])')
        self.upper = eval(pru.search(data).group('ar'))
        # search for Dimensions
        prd = re.compile(r'Dimensions\s*=\s*(?P<rn>\d)')
        self.dims = eval(prd.search(data).group('rn'))
        # search for time range for which sim was run
        prt = re.compile(r'Time\s*=\s*(?P<ar>\[.*\])')
        self.time = eval(prt.search(data).group('ar'))
        # search how many output files were written
        pro = re.compile(r'Out\s*=\s*(?P<rn>\d+)')
        self.out = eval(pro.search(data).group('rn'))

        # open log file and read contents
        try:
            data = open(inp + '_0.log').read()
        except:
            data = open(inp + '.log').read()

        # find how many equations are in the system: this information
        # is impossible to determine from the input file itself and
        # hence we need to parse the log file.
        pme = re.compile(r'Equation system has total (?P<ar>\d+) equations')
        self.meqn = eval(pme.search(data).group('ar'))

        # time between frames
        self.tframe = (self.time[1] - self.time[0])/self.out

        # reduce cell, upper and lower to match dimensions
        self.cells = self.cells[0:self.dims]
        self.lower = self.lower[0:self.dims]
        self.upper = self.upper[0:self.dims]

        # grid spacing
        self.dx = []
        self.axis = []
        for i in range(self.dims):
            dx = (self.upper[i] - self.lower[i])/self.cells[i]
            self.dx.append( dx ) # add to cell size list
            axis = numpy.zeros((self.cells[i],), numpy.float)
            for ix in range(self.cells[i]):
                axis[ix] = self.lower[i] + (ix+0.5)*dx
            self.axis.append( axis ) # add to axis list

        self.extra = extra
        
    def read(self, frame):
        r"""read(frame : int) -> numpy array

        Reads data from the given ``frame`` and returns it as a
        numpy array. The returned data has shape (nx, meqn) for 1D
        simulations, and (nx, ny, meqn) for 2D simulations"""

        # construct file names
        fn = self.base + "_" + self.block + "_" + str(frame) + ".txt"
        print("Reading frame data from file %s...") % fn        
        data = open(fn, 'r').read().split()
        condata = numpy.array( map(lambda s: float(s), data),
                                 numpy.float)
        # now get it into correct shape depending on dimension
        if self.dims == 1:
            # 1D data file
            data = numpy.reshape(condata, (self.cells[0], self.meqn*self.extra))
        elif self.dims == 2:
            # 2D data file
            data = numpy.reshape(condata, (self.cells[0]*self.cells[1], self.meqn*self.extra))
            data = numpy.reshape(data, (self.cells[1], self.cells[0], self.meqn*self.extra))
            data = numpy.transpose(data, (1,0,2))

        return data
