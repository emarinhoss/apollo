#!/usr/bin/env python
r"""Provides a general interface to plot warpx results quickly.

"""

import os
from pylab import *
import wxdata
import wxinteractiveplot

from optparse import OptionParser

# set command line options
parser = OptionParser()
parser.add_option('-i', '--input', action = 'store',
                  dest = 'inputFile',
                  help = 'Base name of simulation')
parser.add_option('-f', '--frame', action = 'store',
                  dest = 'frame',
                  help = 'Frame to plot',
                  default = 0)
parser.add_option('-v', '--variable', action = 'store',
                  dest = 'variable',
                  help = 'Variable to plot',
                  default = 'qnew')
parser.add_option('-c', '--component', action = 'store',
                  dest = 'component',
                  help = 'Component to plot')
parser.add_option('-d', '--dg', action = 'store_true',
                  dest = 'dg',
                  help = 'Plot DG results')
parser.add_option('-m', '--meqn', action='store',
                  dest = 'meqn',
                  help = 'Number of equations')
parser.add_option('--show-variables', action = 'store_true',
                  dest = 'showVariables',
                  help = 'Show variables in output file')
parser.add_option('--variable-info', action = 'store',
                  dest = 'variableInfo',
                  help = 'Show information about variable')
parser.add_option('--save-png', action = 'store_true',
                  dest = 'savePng',
                  help = 'Save png of plot displayed')
parser.add_option('--contour', action = 'store_true',
                  dest = 'contour',
                  help = 'Create a contour plot for 2D data.'
                  ' The default it to make pcolor plots.')
parser.add_option('--dont-show', action = 'store_true',
                  dest = 'dontShow',
                  help = 'Do not show plot',
                  default = False)
parser.add_option('-l', '--plot-limits', action = 'store',
                  dest = 'plotLimits',
                  help = 'Limits for plotting. Set using "l,u" where l is'
                  ' lower limit and u is upper limit')
parser.add_option('-t', '--transform-variable', action = 'store',
                  dest = 'transformVariable',
                  help = 'Name of transform variable to plot')
parser.add_option('--transforms-file', action = 'store',
                  dest = 'transformsFile',
                  help=  'File for variable transforms')
parser.add_option('--show-transform-variables', action='store_true',
                  dest = 'showTransfromVariables',
                  help = 'Show list of transform variables')
parser.add_option('-n', '--interactive', action='store_true',
                  dest = 'interactive',
                  help = 'Start interactive mode')

(options, args) = parser.parse_args()

def make1Dplot(q, options):
    if options.component:
        component = int(options.component)
        data = q[:,component]

        # construct file and ylabel
        fn = '%s_%s_c%d_%.5d.png' % (options.inputFile, options.variable, component, frame)
        yl = '%s[%d]' % (options.variable, component)
    elif options.transformsFile:
        # compile and evaluate transforms file
        sys.path.append(os.path.abspath('.'))
        mod = __import__(options.transformsFile)
        trans = mod.transformregistry[options.transformVariable]
        data = trans(q)

        # construct file and ylabel
        fn = '%s_%s_%.5d.png' % (options.inputFile, options.transformVariable, frame)
        yl = '%s' % options.transformVariable
        pass
    else:
        raise NameError("Must specify component to plot or transform file for variable %s" % options.variable)

    plot(q.grid.X, data)
    title('Time = %g' % d.time)
    ylabel(yl)

    if options.plotLimits:
        lmin, lmax = eval(options.plotLimits)
        gca().set_ylim( (lmin, lmax) )

    # save figure if needed
    if options.savePng:
        savefig(fn)

    # show if needed
    if options.dontShow == False:
        show()

def make2Dplot(q, options):
    if options.component:
        component = int(options.component)
        data = q[:,:,component]

        # construct file and ylabel
        fn = '%s_%s_c%d_%.5d.png' % (options.inputFile, options.variable, component, frame)
        yl = '%s[%d]' % (options.variable, component)
    elif options.transformsFile:
        # compile and evaluate transforms file
        sys.path.append(os.path.abspath('.'))
        mod = __import__(options.transformsFile)
        trans = mod.transformregistry[options.transformVariable]
        data = trans(q)

        # construct file and ylabel
        fn = '%s_%s_%.5d.png' % (options.inputFile, options.transformVariable, frame)
        yl = '%s' % options.transformVariable
        pass
    else:
        raise NameError("Must specify component to plot or transform file for variable %s" % options.variable)
    
    XX, YY = meshgrid(q.grid.X, q.grid.Y)
    if options.contour:
        # make contour plot
        data = data.transpose()
        contour(XX, YY, data)
        axis('image')
    else:
        # make pcolormesh plot
        data = data.transpose()
        if options.plotLimits:
            vmin, vmax = eval(options.plotLimits)
            pcolormesh(XX, YY, data, vmin=vmin, vmax=vmax, shading='flat')
        else:
            pcolormesh(XX, YY, data, shading='flat')

        axis('image')
        hot()
        try:
            colorbar()
        except:
            # this means something is wrong with colorbar
            pass
    title('Time = %g. Variable %s' % (d.time, yl))

    # save figure if needed
    if options.savePng:
        savefig(fn)

    # show if needed
    if options.dontShow == False:
        show()

# basic parameters
frame = int(options.frame)
# open data file to plot
d = wxdata.WxData(options.inputFile, frame)

# print list of variable if requested
if options.showVariables:
    print(d.variables())
    d.close()
    exit()

# print list of transform variables
if options.showTransfromVariables:
    # compile and evaluate transforms file
    sys.path.append(os.path.abspath('.'))
    mod = __import__(options.transformsFile)
    print(mod.transformregistry.keys())
    d.close()
    exit()

# print information about variable if requested
if options.variableInfo:
    q = d.read(options.variableInfo)
    print('Variable "%s" lives on grid "%s" with ') % (options.variableInfo, q.onGrid)
    print(q.grid)
    print('No of components is '), q.numComponents
    d.close()
    exit()

# read variable specified
if options.dg:
    meqn = int(options.meqn)    
    q = d.readDG(options.variable, meqn)
else:
    q = d.read(options.variable)

try:
    if q.grid.ndims == 1:
        make1Dplot(q, options)
    elif q.grid.ndims == 2:
        make2Dplot(q, options)
    elif q.grid.ndims == 3:
        raise Exception("3D plots not currently supported")
except NameError as strerror:
    print(strerror)
if options.interactive:
    ip = wxinteractiveplot.WxInteractivePlot(options)
    ip.cmdloop()

# close HDF5 file
d.close()

