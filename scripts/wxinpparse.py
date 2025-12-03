#!/usr/bin/env python
import os
import wxinputparser
import wxmacroexpander

from optparse import OptionParser

if __name__ == '__main__':
    # read command line options
    parser = OptionParser()
    parser.add_option('-i', '--input', action = 'store',
                      dest = 'inputFile',
                      help = 'Input file to preprocess')
    parser.add_option('-o', '--output', action = 'store',
                      dest = 'outputFile',
                      help = 'Name of preprocessed input file')
    parser.add_option('-x', '--execute', action = 'store',
                      dest = 'program',
                      help = 'Program to execute. '
                      'For example mp will run macroexpansion followed by preprocessor'
                      ' while m will only run the macroexpanion. Defaults to mp.')

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

    # determine name of output file
    if options.outputFile:
        outFile = options.outputFile
    else:
        n = inpFile.split('.')
        outFile = n[0] + '.inp'
        if outFile == inpFile:
            outFile = outFile + '_1'

    # determine sequence of commands to execute
    if options.program:
        program = options.program
    else:
        program = "mp"

    # run command sequence
    tempFileA = inpFile + "_temp1"
    tempFileB = inpFile + "_temp2"

    # open input file and copy it into tempFileA
    d = open(inpFile, "r").read()
    # find the initial <apollo> or <warpx> tag (warpx is legacy name)
    loc = d.find("<apollo>")
    if loc < 0:
        loc = d.find("<warpx>")
    if loc < 0:
        raise Exception("Unable to find opening apollo or warpx tag in input file %s" % inpFile)
    prefix = d[0:loc]

    # add tags to feed to macroexpansion function
    newData = "%s\n#begin python\n%s\n#end python\n%s" % (prefix, prefix, d[loc:])
    # write this data to temporary file
    ta = open(tempFileA, "w")
    ta.write(newData)
    ta.close()
    
    for c in program[:-1]:
        if c=='p':
            print("Running preprocessor")
            wxinputparser.preprocessFile(tempFileA, tempFileB)
        elif c=='m':
            print("Running macroexpansion")
            wxmacroexpander.expandFile(tempFileA, tempFileB)
        else:
            print("Unknown command %s. Skipping ..." % c)

        # switch files
        tf = tempFileA
        tempFileA = tempFileB
        tempFileB = tf

    # final command
    c = program[-1]
    if c=='p':
        print("Running preprocessor")
        wxinputparser.preprocessFile(tempFileA, outFile)
    elif c=='m':
        print("Running macroexpansion")
        wxmacroexpander.expandFile(tempFileA, outFile)
    else:
        print("Unknown command %s. Skipping ..." % c)
