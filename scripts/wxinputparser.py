#!/usr/bin/env python
r"""Input file parser. This code evaluates the top level python code in
each input file and performs the substitution in the rest of the input
file.
"""

import os.path
import re
import xml.dom.minidom
import xml.parsers.expat
from optparse import OptionParser

class WxDataSet(object):
    r"""Stores data in a WarpX input file block
    """
    
    def __init__(self, name):
        r"""
        WxDataSet(name : string) -> WxDataSet
        
        Create a new empty data set with supplied name
        """
        self.name = name # name of dataset
        self.values = {} # name value pairs for names in this dataset
        self.children = [] # list of children datasets
        self.childrenMap = {} # map of name to children dataset
        self.xmlStr = '' # XML string for this data-set
        self.processedValues = [] # key value pairs (processed)
        self.unprocessedValues = [] # key value pairs (unprocesses)

class WxInpParse(object):
    r"""Parses an input file into a fully evaluated WxDataSet object.
    """

    # for removing comments
    replaceCommentsRe = re.compile('#.*')
    # for line continuation
    lineContRe = re.compile(rr"\\\s*")
    # for parsing key = value lines
    nameValueRe = re.compile(rr'\s*(?P<key>[a-zA-Z_][a-zA-Z0-9\-_]*)\s*=\s*(?P<value>.*)')
    # for parsing strings
    strRe = re.compile(r'(?P<str>".*")');

    def __init__(self, inpFile=None):
        if inpFile:
            if not os.path.exists(inpFile):
                raise Exception("File %s does not exist!" % inpFile)
            self.inpFile = inpFile
            self.inp = open(inpFile, 'r')
        self.pyDict = {} # dictionary of names to python compiled objects

        # create top level dataset
        if inpFile:
            self.dataSet = self._parse()

    def _parse(self):
        data = self.inp.read()
        # append a start and end tag before feeding to the XML parser
        data = '<__top>%s</__top>' % data
        # replace comments with empty strings
        ndata, nr = self.replaceCommentsRe.subn('\n', data)
        # flatten lines with continuation symbol '\' at the end
        ndata, nr = self.lineContRe.subn('', ndata)
        # add additional quotes around all strings
        def addQuotes(mo):
            newStr = "'%s'" % mo.group('str')
            return newStr
        ndata, nr = self.strRe.subn(addQuotes, ndata)
        
        try:
            topNode = xml.dom.minidom.parseString(ndata).childNodes[0]
        except xml.parsers.expat.ExpatError as e:
            # print error
            lineno = e.lineno
            lines = ndata.split('\n')
            closeTag = lines[lineno-1]
            print("** Error: Closing tag %s does not match any in the input file" % (closeTag.strip()))
            exit(1)
                
        code = '' # python code string in topNode
        elemNodes = [] # element nodes in topNode
        # parse the text nodes of top node into a string
        for child in topNode.childNodes:
            if child.nodeType == child.TEXT_NODE:
                code = code + child.wholeText
            elif child.nodeType == child.ELEMENT_NODE:
                elemNodes.append(child)

        # compile and evaluate the code
        co = compile(code, self.inpFile + '.err', 'exec')
        exec(co)
        # add all symbols in co to self
        for name in co.co_names:
            try:
                self.pyDict[name] = eval(name)
            except:
                # do nothing if adding failed
                pass
        
        # recursively step into the XML file creating WarpX data object
        wxds = self._createDataSet(elemNodes[0]) # there is only one top level block

        # create map of name -> children
        self._makeChildrenMap(wxds)
        
        return wxds

    def _makeChildrenMap(self, dataset):
        for c in dataset.children:
            dataset.childrenMap[c.name] = c
            self._makeChildrenMap(c)

    def _createDataSet(self, domElem):
        wxds = WxDataSet(domElem.tagName)

        # walk over node's text, and store it for substitution
        code = ''
        for child in domElem.childNodes:
            if child.nodeType == child.TEXT_NODE:
                code = code + child.wholeText

        # split code into name-value pairs
        moItr = self.nameValueRe.finditer(code)
        for mo in moItr:
            key = mo.group('key')
            valueStr = mo.group('value')

            # add unprocessed pair to dataset
            wxds.unprocessedValues.append( (key, valueStr) )

            # now parse value into tokens, performing needed substitution
            co = compile(valueStr, self.inpFile + '.err', 'exec')
            # create dictionary of names for default replacement
            di = {}
            for name in co.co_names:
                di[name] = name
            # now perform replacement and add it to the dataset
            try:
                wxds.values[key] = eval(valueStr, di, self.pyDict)
                # add processed pair to dataset
                wxds.processedValues.append( (key, wxds.values[key]) )
            except:
                # if evaluation failed, assume the rhs is a string
                wxds.values[key] = valueStr
                wxds.processedValues.append( (key, valueStr) )
            
        # walk over node's children, adding them to children list
        for child in domElem.childNodes:
            if child.nodeType == child.ELEMENT_NODE:
                child_ds = self._createDataSet(child)
                wxds.children.append(child_ds)

        return wxds

    def toProcessedXml(self):
        res = self._dataSetToXml(self.dataSet)
        repl = re.compile("'")
        newRes, nr = repl.subn('', res)
        return newRes

    def _dataSetToXml(self, dataset, spaces=2):
        # for each name=value pair insert string into XML
        myValues = ''
        for pair in dataset.processedValues:
            myValues = '%s\n%s%s = %s' % (myValues, spaces*' ', pair[0], str(pair[1]))

        xmlString = ''
        # now get each child's Xml and append it to ours
        for child in dataset.children:
            xmlString = '%s\n%s' % (xmlString,
                                    self._dataSetToXml(child, spaces = spaces + 2))

        # now attach start and end tags and return
        spt = (spaces-2)*' '
        return '%s<%s>%s\n%s\n%s</%s>\n' % (spt, dataset.name, myValues, xmlString, spt, dataset.name)

def preprocessFile(inFileName, outFileName):
    wip = WxInpParse(inFileName)
    d = wip.toProcessedXml()
    fp = open(outFileName, "w")
    fp.writelines(d)
    fp.close()
