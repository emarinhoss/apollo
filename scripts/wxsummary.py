#!/usr/bin/env python
import os
import sys
import wxinputparser

scriptdir = os.path.split( os.path.realpath( sys.argv[0] ) )[0]

from optparse import OptionParser

# read command line options
parser = OptionParser()
parser.add_option('-i', '--input', action = 'store',
                  dest = 'inputFile',
                  help = 'Input file to preprocess')
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

n = inpFile.split('.')
outFile = n[0] + '.html'

# run preprocessor to create inp file
p = os.popen('%s/wxinpparse.py -i %s' % (scriptdir, inpFile))
p.read()

# open file for reading
pwd = fn = os.path.abspath('.')
inp = wxinputparser.WxInpParse('%s/%s.inp' % (pwd, n[0]))
# get top level object
wx = inp.dataSet.children[0]
# create directory for storing HTML for each subsolver
try:
    os.mkdir('%s__data' % n[0])
except:
    # directory exists so just keep going
    pass

# now loop over input file elements and create one subsolver file per
# subsolver
for elem in wx.childrenMap:
    e = wx.childrenMap[elem]
    if 'Type' in e.values:
        if e.values['Type'] == 'WxSubSolver':
            txt = inp._dataSetToXml(e) # text of subsolver
            # replace < and > characters with proper HTML symbols
            txt = txt.replace('<', '&lt;')
            txt = txt.replace('>', '&gt;')

            # open a new file and write data out            
            fp = open('%s__data/%s.html' % (n[0], elem), 'w')
            fp.writelines('<html><body>\n<pre>\n%s\n</pre></body></html>' % txt)
            fp.close()

# get subsolver sequence
ss = wx.childrenMap['SolverSequence']

# open HTML file for output
hf = open(outFile, 'w')
# write header information, including CSS
header = r"""
<?xml version="1.0" encoding="utf-8" ?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
<head>
<SCRIPT TYPE="text/javascript">
<!--
function popup(mylink, windowname)
{
if (! window.focus)return true;
var href;
if (typeof(mylink) == 'string')
   href=mylink;
else
   href=mylink.href;
window.open(href, windowname, 'width=400,height=200,scrollbars=yes');
return false;
}
//-->
</SCRIPT>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<meta name="generator" content="Docutils 0.5: http://docutils.sourceforge.net/" />
<title>Summary of input file %s</title>
<style type="text/css">
%s
</style>
</head>
<body>
""" % (n[0], open('%s/default.css' % scriptdir).read())
hf.writelines(header)
hf.writelines("""<h1 class="title">Summary of input file %s</h1>""" % inpFile)

# write list of subsolver steps
hf.writelines("""<h2 class="title">StartOnly subsolver steps</h2>""")

hf.writelines("<p>")
hf.writelines("""<table class="mctable" cellpadding="3" cellspacing="0" >""")
hf.writelines("<tr>")
hf.writelines("""<td class="mctablehead">Step</td>""")
hf.writelines("""<td class="mctablehead">SubSolvers</td>""")
hf.writelines("""<td class="mctablehead">DtFrac</td>""")
hf.writelines("""<td class="mctablehead">SyncVars</td>""")
hf.writelines("</tr>")

# list of subsolvers
completeSubSolverLst = []
if 'StartOnly' in ss.values:
    startOnly = ss.values['StartOnly']
else:
    startOnly = []
# loop over each step and extract subsolvers in that step
for step in startOnly:
    # get solver step
    solverStep = wx.childrenMap[step]
    # from this extract names of subsolvers
    subSolvers = solverStep.values['SubSolvers']
    completeSubSolverLst.extend(subSolvers)
    # time-fraction
    if 'DtFrac' in solverStep.values:
        dt = solverStep.values['DtFrac']
    else:
        dt = '1.0'
    # SyncVars
    if 'SyncVars' in solverStep.values:
        syncVars = solverStep.values['SyncVars']
    else:
        syncVars = ['-']

    ssString = ''
    for ss in subSolvers:
        ns = '<a href="%s__data/%s.html" onClick="return popup(this, \'%s\')">%s</a>' % (n[0], ss, ss, ss)
        ssString = '%s %s' % (ssString, ns)

    hf.writelines("<tr>")
    hf.writelines("<td>%s</td>" % step)
    hf.writelines("<td>%s</td>" % ssString)
    hf.writelines("<td>%s</td>" % dt)
    hf.writelines("<td>%s</td>" % ', '.join(syncVars))

    hf.writelines("</tr>")
hf.writelines('</table>')

# write details for each subsolver
hf.writelines("""<h2 class="title">StartOnly subsolver details</h2>""")

hf.writelines("<p>")
hf.writelines("""<table class="mctable" cellpadding="3" cellspacing="0" >""")
hf.writelines("<tr>")
hf.writelines("""<td class="mctablehead">SubSolver</td>""")
hf.writelines("""<td class="mctablehead">Kind</td>""")
hf.writelines("""<td class="mctablehead">ReadVars</td>""")
hf.writelines("""<td class="mctablehead">WriteVars</td>""")
hf.writelines("</tr>")

for subsolver in completeSubSolverLst:
    # get hold of object
    ss = wx.childrenMap[subsolver]
    # get kind field
    if 'Kind' in ss.values:
        kind = ss.values['Kind']
    else:
        kind = ['-']    

    # get read vars
    if 'ReadVars' in ss.values:
        readVars = ss.values['ReadVars']
    else:
        readVars = ['-']

    # get write vars
    if 'WriteVars' in ss.values:
        writeVars = ss.values['WriteVars']
    else:
        writeVars = ['-']

    ns = '<a href="%s__data/%s.html" onClick="return popup(this, \'%s\')">%s</a>' \
         % (n[0], subsolver, subsolver, subsolver)
    hf.writelines("<tr>")
    hf.writelines("<td>%s</td>" % ns)
    hf.writelines("<td>%s</td>" % kind)
    hf.writelines("<td>%s</td>" % ', '.join(readVars))
    hf.writelines("<td>%s</td>" % ', '.join(writeVars))

    hf.writelines("</tr>")        
hf.writelines('</table>')

# write list of subsolver steps
hf.writelines("""<h2 class="title">PerStep subsolver steps</h2>""")

hf.writelines("<p>")
hf.writelines("""<table class="mctable" cellpadding="3" cellspacing="0" >""")
hf.writelines("<tr>")
hf.writelines("""<td class="mctablehead">Step</td>""")
hf.writelines("""<td class="mctablehead">SubSolvers</td>""")
hf.writelines("""<td class="mctablehead">DtFrac</td>""")
hf.writelines("""<td class="mctablehead">SyncVars</td>""")
hf.writelines("</tr>")

# list of subsolvers
ss = wx.childrenMap['SolverSequence']
completeSubSolverLst = []
if 'PerStep' in ss.values:
    perStep = ss.values['PerStep']
else:
    perStep = []
# loop over each step and extract subsolvers in that step
for step in perStep:
    # get solver step
    solverStep = wx.childrenMap[step]
    # from this extract names of subsolvers
    subSolvers = solverStep.values['SubSolvers']
    completeSubSolverLst.extend(subSolvers)
    # time-fraction
    if 'DtFrac' in solverStep.values:
        dt = solverStep.values['DtFrac']
    else:
        dt = '1.0'
    # SyncVars
    if 'SyncVars' in solverStep.values:
        syncVars = solverStep.values['SyncVars']
    else:
        syncVars = ['-']

    ssString = ''
    for ss in subSolvers:
        ns = '<a href="%s__data/%s.html" onClick="return popup(this, \'%s\')">%s</a>' % (n[0], ss, ss, ss)
        ssString = '%s %s' % (ssString, ns)

    hf.writelines("<tr>")
    hf.writelines("<td>%s</td>" % step)
    hf.writelines("<td>%s</td>" % ssString)
    hf.writelines("<td>%s</td>" % dt)
    hf.writelines("<td>%s</td>" % ', '.join(syncVars))

    hf.writelines("</tr>")
hf.writelines('</table>')

# write details for each subsolver
hf.writelines("""<h2 class="title">PerStep subsolver details</h2>""")

hf.writelines("<p>")
hf.writelines("""<table class="mctable" cellpadding="3" cellspacing="0" >""")
hf.writelines("<tr>")
hf.writelines("""<td class="mctablehead">SubSolver</td>""")
hf.writelines("""<td class="mctablehead">Kind</td>""")
hf.writelines("""<td class="mctablehead">ReadVars</td>""")
hf.writelines("""<td class="mctablehead">WriteVars</td>""")
hf.writelines("</tr>")

for subsolver in completeSubSolverLst:
    # get hold of object
    ss = wx.childrenMap[subsolver]
    # get kind field
    if 'Kind' in ss.values:
        kind = ss.values['Kind']
    else:
        kind = ['-']    

    # get read vars
    if 'ReadVars' in ss.values:
        readVars = ss.values['ReadVars']
    else:
        readVars = ['-']

    # get write vars
    if 'WriteVars' in ss.values:
        writeVars = ss.values['WriteVars']
    else:
        writeVars = ['-']

    ns = '<a href="%s__data/%s.html" onClick="return popup(this, \'%s\')">%s</a>' \
         % (n[0], subsolver, subsolver, subsolver)
    hf.writelines("<tr>")
    hf.writelines("<td>%s</td>" % ns)
    hf.writelines("<td>%s</td>" % kind)
    hf.writelines("<td>%s</td>" % ', '.join(readVars))
    hf.writelines("<td>%s</td>" % ', '.join(writeVars))

    hf.writelines("</tr>")        
hf.writelines('</table>')

hf.writelines('</body></html>')
