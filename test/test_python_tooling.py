#!/usr/bin/env python3
"""Regression tests for the Python tooling under scripts/.

These exist because an automated Python 2-to-3 pass once left every script in
the directory either unparseable or guaranteed to raise on its first print.
Each test below pins one of the failure modes that pass introduced, so a repeat
is caught before it is committed rather than the next time somebody tries to
post-process a run.

Runs with the standard library alone; pyflakes adds one more check when present.

    python3 -m unittest discover -s test -v
"""

import ast
import io
import os
import re
import subprocess
import sys
import tokenize
import unittest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = os.path.join(REPO, 'scripts')
BUILDCONF = os.path.join(REPO, 'src', 'buildconf')


def python_files():
    """Every Python file this repository is responsible for."""
    found = []
    for root in (SCRIPTS, BUILDCONF, os.path.join(REPO, 'test')):
        for dirpath, dirnames, filenames in os.walk(root):
            dirnames[:] = [d for d in dirnames if d != '__pycache__']
            found += [os.path.join(dirpath, f)
                      for f in filenames if f.endswith('.py')]
    found.append(os.path.join(REPO, 'src', 'this_host_config.py'))
    return sorted(p for p in found if os.path.isfile(p))


def rel(path):
    return os.path.relpath(path, REPO)


def code_lines(path):
    """Yield (lineno, text) for lines that are not comments or string bodies.

    Tokenising rather than grepping keeps the checks from firing on commented-out
    legacy code or on text that merely appears inside a docstring.
    """
    with open(path, 'rb') as fh:
        source = fh.read()
    try:
        tokens = list(tokenize.tokenize(io.BytesIO(source).readline))
    except (tokenize.TokenError, SyntaxError, IndentationError):
        return
    skip = set()
    for tok in tokens:
        if tok.type in (tokenize.COMMENT, tokenize.STRING):
            for n in range(tok.start[0], tok.end[0] + 1):
                skip.add(n)
    for lineno, text in enumerate(source.decode('utf-8').splitlines(), 1):
        if lineno not in skip:
            yield lineno, text


class TestEveryFileParses(unittest.TestCase):
    """scripts/SGC.py and three others once failed at import with TabError."""

    def test_parses(self):
        broken = []
        for path in python_files():
            with open(path, encoding='utf-8') as fh:
                src = fh.read()
            try:
                ast.parse(src, filename=path)
            except SyntaxError as exc:
                broken.append(f'{rel(path)}:{exc.lineno}: {exc.msg}')
        self.assertEqual([], broken, '\n'.join(['files that do not parse:'] + broken))


class TestNoMangledPrints(unittest.TestCase):
    """`print(A) + B` parses but raises TypeError the moment it runs.

    Python 2's `print A + B` was mechanically rewritten to `print(A) + B`, which
    tries to add a str to the None that print() returns.
    """

    PATTERN = re.compile(r'^\s*print\s*\(.*\)\s*\+')

    def test_no_dangling_concatenation(self):
        offenders = []
        for path in python_files():
            for lineno, text in code_lines(path):
                if self.PATTERN.match(text):
                    offenders.append(f'{rel(path)}:{lineno}: {text.strip()}')
        self.assertEqual([], offenders, '\n'.join(
            ['`print(A) + B` raises TypeError at runtime; write `print(A + B)`:']
            + offenders))


class TestNoRemovedNumpyAliases(unittest.TestCase):
    """numpy.float and friends were removed in NumPy 1.24 (December 2022)."""

    PATTERN = re.compile(
        r'\b(?:numpy|np)\.(float|int|bool|object|str|complex|long|unicode)\b'
        r'(?![\w.])')

    REPLACEMENT = {
        'float': 'float64 (or the builtin float)',
        'int': 'int64 (or the builtin int)',
        'bool': 'bool_ (or the builtin bool)',
        'object': 'object_ (or the builtin object)',
        'str': 'str_ (or the builtin str)',
        'complex': 'complex128 (or the builtin complex)',
        'long': 'int64',
        'unicode': 'str_',
    }

    def test_no_removed_aliases(self):
        offenders = []
        for path in python_files():
            for lineno, text in code_lines(path):
                m = self.PATTERN.search(text)
                if m:
                    offenders.append(
                        f'{rel(path)}:{lineno}: numpy.{m.group(1)} was removed in '
                        f'NumPy 1.24; use numpy.{self.REPLACEMENT[m.group(1)]}')
        self.assertEqual([], offenders, '\n'.join(offenders))


class TestNoPython2Builtins(unittest.TestCase):
    """Names that exist only in Python 2 and raise NameError under Python 3."""

    NAMES = ('raw_input', 'xrange', 'unicode', 'basestring', 'execfile', 'reduce')

    def test_no_python2_names(self):
        offenders = []
        for path in python_files():
            with open(path, encoding='utf-8') as fh:
                try:
                    tree = ast.parse(fh.read(), filename=path)
                except SyntaxError:
                    continue  # reported by TestEveryFileParses
            imported = set()
            for node in ast.walk(tree):
                if isinstance(node, ast.ImportFrom):
                    imported.update(a.asname or a.name for a in node.names)
                elif isinstance(node, ast.Import):
                    imported.update((a.asname or a.name).split('.')[0]
                                    for a in node.names)
            for node in ast.walk(tree):
                if (isinstance(node, ast.Name) and isinstance(node.ctx, ast.Load)
                        and node.id in self.NAMES and node.id not in imported):
                    offenders.append(
                        f'{rel(path)}:{node.lineno}: {node.id} does not exist in '
                        f'Python 3')
        self.assertEqual([], offenders, '\n'.join(offenders))


class TestNoMixedIndentation(unittest.TestCase):
    """A file that indents with both tabs and spaces raises TabError.

    Python 2 resolved a tab to the next 8-column stop and accepted the mix;
    Python 3 refuses the file outright, so this is an import-time failure rather
    than a style nit. Indenting consistently with tabs is fine, so the test looks
    at the INDENT tokens the tokenizer actually produces instead of at raw
    leading whitespace, which would also flag continuation lines.
    """

    def test_no_mixed_indentation(self):
        offenders = []
        for path in python_files():
            with open(path, 'rb') as fh:
                data = fh.read()
            try:
                tokens = list(tokenize.tokenize(io.BytesIO(data).readline))
            except (tokenize.TokenError, SyntaxError, IndentationError):
                continue  # reported by TestEveryFileParses
            indents = [t.string for t in tokens if t.type == tokenize.INDENT]
            tabbed = any('\t' in i for i in indents)
            spaced = any(' ' in i for i in indents)
            if tabbed and spaced:
                first = next(t for t in tokens
                             if t.type == tokenize.INDENT and '\t' in t.string)
                offenders.append(
                    f'{rel(path)}: indents with both tabs and spaces '
                    f'(first tab indent at line {first.start[0]})')
        self.assertEqual([], offenders, '\n'.join(
            ['mixing tab and space indentation raises TabError:'] + offenders))


class TestNoUndefinedNames(unittest.TestCase):
    """pyflakes catches typo'd identifiers, the class of bug that hid
    WxDGArray2, Console() and open(fileName) in otherwise valid files."""

    # src/buildconf/*.py are SConscript fragments, not standalone modules: SCons
    # execs them with these names already bound, so pyflakes cannot see them.
    SCONS_GLOBALS = ('warpMConstructionEnv', 'Import', 'Export', 'Configure',
                     'Exit', 'SConscript', 'Environment', 'Return', 'Help',
                     'Variables', 'BoolVariable', 'SetOption', 'GetOption',
                     'env')

    def test_pyflakes_clean(self):
        try:
            import pyflakes  # noqa: F401
        except ImportError:
            self.skipTest('pyflakes not installed (pip install pyflakes)')
        result = subprocess.run(
            [sys.executable, '-m', 'pyflakes'] + python_files(),
            capture_output=True, text=True)
        sconsNames = {f"undefined name '{n}'" for n in self.SCONS_GLOBALS}
        hard = []
        for line in result.stdout.splitlines():
            if 'unable to detect undefined names' in line:
                # `from numpy import *` leaves pyflakes unable to resolve names at
                # all; that advisory is not itself an undefined-name finding.
                continue
            if 'undefined name' not in line:
                continue
            if 'may be undefined, or defined from star imports' in line:
                continue
            if any(n in line for n in sconsNames):
                continue
            hard.append(line.replace(REPO + os.sep, ''))
        self.assertEqual([], hard, '\n'.join(['undefined names:'] + hard))


if __name__ == '__main__':
    unittest.main()
