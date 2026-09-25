#!/usr/bin/env python
import sys
import os
import time
from optparse import OptionParser

start_time = time.perf_counter()

class PymException(Exception): pass
class PymEndOfFile(PymException): pass
class PymExit(PymException): pass

PYM_PATH = []
PYM_PREFIX_MAP = { '/': ("/","END_") }  # This translates a prefix so that
    # special macro names can be used:
    # e.g <[/KEY]> => <[END_KEY]>
    # Prefixes must be recognizable by
    # their 1st character (key in map).
PYM_EXPRESSION = ["<[","]>"]

ENVIRONMENT = {
   "PYM_EXPRESSION": PYM_EXPRESSION,
   "PYM_PREFIX_MAP": PYM_PREFIX_MAP,
   "PYM_PATH": PYM_PATH,
   "PymException": PymException,       ## mirror the exceptions in pym
   "PymEndOfFile": PymEndOfFile,       ## environment
   "PymExit": PymExit,
}

HTTP_HEADER = """Content-type: text/html

"""

HTTP_FOOTER = """
<!-- pym timings: %g sec rendering, %g sec overall -->
"""

def pym_determine_pym_hash(line):
    # check if the line is special line or just a python comment line
    if line.find("end python") > 0:
        return True
    if line.find("begin python") > 0:
        return True
    if line.find("include") == 1:
        return True
    if line.find("include_direct") == 1:
        return True
    if line.find("if") == 1:
        return True
    if line.find("elif") == 1:
        return True
    if line.find("else") == 1:
        return True
    if line.find("endif") == 1:
        return True

    # just a python comment
    return False

def pym_error(message, loc):
    print("ERROR:", message, "in '%s'.", loc[0])
    sys.exit(-1)

def pym_expand(text, env, loc, out):
    (begin,end) = PYM_EXPRESSION
    prefix_map = PYM_PREFIX_MAP
    begin_len = len(begin)
    end_len = len(end)
    if len(text) <= begin_len + end_len:
        out.append(text)
    else:
        pos = 0
        start = text.find(begin, pos)
        while (start >= 0):
            stop = text.find(end, start+begin_len)
            if stop < 0: pym_error("unterminated python macro", loc)
            out.append(text[pos:start])
            exp = text[start+begin_len:stop].strip()
            prefix = prefix_map.get(exp[0])
            try:
                if prefix:
                    len_pr = len(prefix[0])
                    if len(exp) >= len_pr and exp[:len_pr] == prefix[0]:
                        value = eval(prefix[1]+exp[len_pr:], env, env)
                    else:
                        pym_error("illegal prefix '%s'" % exp[:len_pr], loc)
                else:
                    value = eval(exp, env, env)
            except PymException: raise
            except NameError: raise     # doesn't have lineno
            except KeyError: raise      # ditto
            except ImportError: raise           # ..
            except AttributeError: raise        # ditto
            except TypeError: raise     # ditto
            except Exception as error:
                error.filename = loc[0]
                error.lineno = error.lineno + loc[1] + \
                               len(text[0:start].split('\n'))
                raise
            pym_expand(str(value), env, loc, out)
            pos = stop+end_len
            start = text.find(begin, pos)
        out.append(text[pos:])

def pym_read_file(filename, out):
    file = open(filename,"r")
    out.append(file.read())
    file.close()

def pym_expand_file(filename, env, out):
    file = open(filename,"r")
    text = file.read()
    file.close()
    lnum = 1
    py_pos = -1
    tx_pos = 0
    loc = (filename, 0)
    pos = 0
    cond = 1
    condstack = []
    lines = text.split('\n')
    if len(lines[0]) > 2 and lines[0][:2] == "#!":
        lines = lines[1:]
        lnum = 2
    for line in lines:
        end = pos + len(line) + 1
        if line and line[0] == '#':
            if pym_determine_pym_hash(line):
                if tx_pos >= 0:
                    tx_start = tx_pos ; tx_pos = -1
                    try:
                        if cond:
                            pym_expand(text[tx_start:pos], env, loc, out)
                    except PymEndOfFile:
                        break
                if line.find("end python") > 0:
                    if py_pos < 0: pym_error("superfluous end python", loc)
                    if cond:
                        try:
                            exec(text[py_pos:pos], env, env)
                        except PymExit: raise
                        except PymEndOfFile:
                            py_pos = -1
                            break
                        except NameError: raise     # doesn't have lineno
                        except KeyError: raise      # ditto
                        except ImportError: raise   # ditto
                        except AttributeError: raise    # ditto
                        except TypeError: raise     # ditto
                        except Exception as error:
                            error.filename = loc[0]
                            error.lineno = error.lineno + loc[1]
                            raise
                    py_pos = -1
                    tx_pos = end
                elif line.find("begin python") > 0:
                    py_pos = end
                    loc = (filename, lnum)
                elif line.find("include") == 1:
                    namestart = 8
                    if line.find("include_direct") == 1:
                        namestart = 15
                    if cond:
                        dir = os.path.dirname(filename)
                        include = eval(line[namestart:].strip(), env, env)
                        includefilename = os.path.join(dir,include)
                        if not os.path.isfile(includefilename):
                            for dir in PYM_PATH:
                                includefilename = os.path.join(dir, include)
                                if os.path.isfile(includefilename): break
                        if namestart == 8:
                            pym_expand_file(includefilename, env, out)
                        else:
                            pym_read_file(includefilename, out)
                    tx_pos = end
                    loc = (filename, lnum)
                elif line.find("if") == 1:
                    condstack.append(cond)
                    cond = eval(line[3:].strip(), env, env)
                    tx_pos = end
                    loc = (filename, lnum)
                elif line.find("elif") == 1:
                    cond = eval(line[5:].strip(), env, env)
                    tx_pos = end
                    loc = (filename, lnum)
                elif line.find("else") == 1:
                    cond = not cond
                    tx_pos = end
                    loc = (filename, lnum)
                elif line.find("endif") == 1:
                    cond = condstack.pop()
                    tx_pos = end
                    loc = (filename, lnum)
        pos = end
        lnum = lnum + 1

    if py_pos >= 0: pym_error("unterminated python code", loc)
    len_text = len(text)
    if tx_pos >= 0 and tx_pos < len_text:
        try: pym_expand(text[tx_pos:min(pos,len_text)], env, loc, out)
        except PymEndOfFile: pass

def expandFile(inFileName, outFileName):
    env = ENVIRONMENT.copy()
    out = []    
    pym_expand_file(inFileName, env, out)
    data = ''
    for text in out: data = data + text

    outFile = open(outFileName, "w")
    if data.strip() != '':
        # write data to output file
        outFile.write(data)
    else:
        # just copy input file to output file
        data = open(inFileName).read()
        outFile.write(data)

    outFile.close()
