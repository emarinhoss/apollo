#!/usr/bin/env python
r"""Creates a shell script which can scripted to provide easy access to
WarpX variables. In particular, the following environment variables
are defined.

SERWARPX : Serial WarpX executable.
PARWARPX : Parallel WarpX executable.
PREWARPX : WarpX pre-processor.
PLOTWARPX: WarpX plotting script.

"""
import os
import sys

scriptDir = sys.path[0]
warpxTopDir = scriptDir[0:scriptDir.find('scripts')]

pf = open('wxpaths.sh', 'w')
pf.writelines('SERWARPX=%s/src/build/xwarpx/xwarpx\n' % warpxTopDir)
pf.writelines('PARWARPX=%s/src/build-par/xwarpx/xwarpx\n' % warpxTopDir)
pf.writelines('PREWARPX=%s/scripts/wxinpparse.py\n' % warpxTopDir)
pf.writelines('PLOTWARPX=%s/scripts/wxplot.py\n' % warpxTopDir)

pf.close()

