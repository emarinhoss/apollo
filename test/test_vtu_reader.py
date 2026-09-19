#!/usr/bin/env python3
"""Does the VTU reader handle the layouts PETSc actually writes?

Apollo's .vtu files are written by PETSc, so their layout follows the PETSc
version rather than anything in this repository:

    PETSc 3.15   no header_type attribute   4-byte length prefix per block
    PETSc 3.19   header_type="UInt64"       8-byte length prefix per block

test/vtu.py assumed the 8-byte form unconditionally. That is correct on Ubuntu
24.04 and wrong on 22.04 - and therefore in Google Colab, where reading a real
run's output died with numpy's "buffer size must be a multiple of element
size", naming neither the file nor the reason.

No PETSc 3.15 is installed here, so these tests BUILD both layouts from the
same arrays and require the reader to return those arrays either way. The
builder is deliberately independent of the reader: it writes the offsets and
prefixes itself, from the format definition, so agreement between them is
evidence rather than tautology.
"""

import os
import struct
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

try:
    import numpy as np
except ImportError:
    np = None

import vtu

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

_PREFIX = {'UInt32': '<I', 'UInt64': '<Q'}


def build_vtu(arrays, header_type='UInt64', declare=True, byte_order='LittleEndian',
              pieces=None):
    """Write a minimal UnstructuredGrid .vtu and return the bytes.

    `arrays` is a list of (name, numpy array). A 2-D array is written with
    NumberOfComponents equal to its second dimension, which is how PETSc writes
    coordinates. `declare=False` omits the header_type attribute, which is what
    PETSc did before 3.16 and what makes UInt32 the operative default.
    """
    fmt = _PREFIX[header_type] if byte_order == 'LittleEndian' \
        else _PREFIX[header_type].replace('<', '>')
    width = struct.calcsize(fmt)

    # `pieces` is a list of array-lists, one per <Piece> - which is how PETSc
    # writes a multi-rank run. The single-piece case is the common one.
    if pieces is None:
        pieces = [arrays]

    piece_xml, blocks, offset = [], [], 0
    for piece_arrays in pieces:
        declarations = []
        npoints = ncells = 1
        for name, data in piece_arrays:
            ncomp = data.shape[1] if data.ndim > 1 else 1
            vtk_type = {'float64': 'Float64', 'float32': 'Float32',
                        'int32': 'Int32', 'int64': 'Int64',
                        'uint8': 'UInt8'}[data.dtype.name]
            declarations.append(
                '        <DataArray type="%s" Name="%s" NumberOfComponents="%d" '
                'format="appended" offset="%d" />'
                % (vtk_type, name, ncomp, offset))
            # The declared byte order applies to the data as well as the prefix,
            # so write it that way rather than leaving native bytes under a
            # BigEndian declaration - that would test the builder, not the reader.
            raw = data.astype(data.dtype.newbyteorder(
                '>' if byte_order == 'BigEndian' else '<')).tobytes()
            blocks.append(struct.pack(fmt, len(raw)) + raw)
            offset += width + len(raw)
            if name == 'Position':
                npoints = len(data)
            elif name == 'types':
                ncells = len(data)
        piece_xml.append(
            '    <Piece NumberOfPoints="%d" NumberOfCells="%d">\n%s\n    </Piece>'
            % (npoints, ncells, '\n'.join(declarations)))

    attrs = 'type="UnstructuredGrid" version="0.1" byte_order="%s"' % byte_order
    if declare:
        attrs += ' header_type="%s"' % header_type

    header = (
        '<?xml version="1.0"?>\n'
        '<VTKFile %s>\n'
        '  <UnstructuredGrid>\n'
        '%s\n'
        '  </UnstructuredGrid>\n'
        '  <AppendedData encoding="raw">\n_'
        % (attrs, '\n'.join(piece_xml)))
    return header.encode('ascii') + b''.join(blocks) + \
        b'\n  </AppendedData>\n</VTKFile>\n'


def write_temp(payload):
    handle = tempfile.NamedTemporaryFile(suffix='.vtu', delete=False)
    handle.write(payload)
    handle.close()
    return handle.name


@unittest.skipIf(np is None, 'numpy not installed')
class TestVtuReader(unittest.TestCase):

    def setUp(self):
        self.arrays = [
            ('Position', np.arange(9, dtype='float64').reshape(3, 3) * 1.5),
            ('connectivity', np.array([0, 1, 2], dtype='int32')),
            ('types', np.array([5], dtype='uint8')),
            ('rhoe', np.array([1.25, -3.5, 7.0], dtype='float64')),
        ]
        self.paths = []

    def tearDown(self):
        for path in self.paths:
            try:
                os.unlink(path)
            except OSError:
                pass

    def _read(self, **kwargs):
        path = write_temp(build_vtu(self.arrays, **kwargs))
        self.paths.append(path)
        return vtu.read(path)

    def _assert_matches(self, got):
        for name, expected in self.arrays:
            self.assertIn(name, got)
            np.testing.assert_array_equal(got[name], expected,
                                          err_msg='array %r differs' % name)

    def test_petsc_319_layout(self):
        """8-byte prefix, header_type declared. What 24.04 ships."""
        self._assert_matches(self._read(header_type='UInt64', declare=True))

    def test_petsc_315_layout(self):
        """4-byte prefix, no header_type attribute. What 22.04 and Colab ship.

        This is the case that failed. The attribute being absent is the whole
        point: the reader has to know the VTK default is UInt32, not assume the
        width it last saw.
        """
        self._assert_matches(self._read(header_type='UInt32', declare=False))

    def test_uint32_declared_explicitly(self):
        """Some writers do declare UInt32. Same 4-byte layout."""
        self._assert_matches(self._read(header_type='UInt32', declare=True))

    def test_big_endian(self):
        """byte_order is read too, not assumed little."""
        self._assert_matches(self._read(header_type='UInt64', declare=True,
                                        byte_order='BigEndian'))

    def test_wrong_width_is_reported_clearly(self):
        """A file whose prefix does not match its declaration must not slip through.

        Before the fix this surfaced as numpy's "buffer size must be a multiple
        of element size", which names neither the file nor the array nor the
        cause. Build the mismatch deliberately: 4-byte prefixes in the data, but
        the header claiming UInt64.
        """
        payload = build_vtu(self.arrays, header_type='UInt32', declare=False)
        payload = payload.replace(b'byte_order="LittleEndian"',
                                  b'byte_order="LittleEndian" header_type="UInt64"', 1)
        path = write_temp(payload)
        self.paths.append(path)
        with self.assertRaises(ValueError) as caught:
            vtu.read(path)
        message = str(caught.exception)
        self.assertIn('header_type', message)
        self.assertIn(os.path.basename(path), message)

    def test_compressed_is_refused_not_guessed(self):
        """Compressed blocks have a different header shape entirely."""
        payload = build_vtu(self.arrays).replace(
            b'byte_order="LittleEndian"',
            b'byte_order="LittleEndian" compressor="vtkZLibDataCompressor"', 1)
        path = write_temp(payload)
        self.paths.append(path)
        with self.assertRaises(ValueError) as caught:
            vtu.read(path)
        self.assertIn('compress', str(caught.exception).lower())

    def test_unknown_header_type_is_refused(self):
        payload = build_vtu(self.arrays).replace(
            b'header_type="UInt64"', b'header_type="UInt16"', 1)
        path = write_temp(payload)
        self.paths.append(path)
        with self.assertRaises(ValueError) as caught:
            vtu.read(path)
        self.assertIn('UInt16', str(caught.exception))

    def test_every_piece_is_read_not_just_the_last(self):
        """PETSc writes one <Piece> per MPI rank; all of them are the file.

        Every Piece re-declares the same array names at its own offset, so a
        reader that walks DataArray elements into a flat dict keeps only the
        last rank's copy. That is what produced the repository's belief that
        "Apollo writes roughly 1/N of the cells on N ranks": the numbers behind
        it were the size of the last piece, not of the file.
        """
        pieces = [
            [('Position', np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]])),
             ('connectivity', np.array([0, 1, 2], dtype='int32')),
             ('offsets', np.array([3], dtype='int32')),
             ('types', np.array([5], dtype='uint8')),
             ('rhoe', np.array([10.0]))],
            [('Position', np.array([[2., 0., 0.], [3., 0., 0.], [2., 1., 0.]])),
             ('connectivity', np.array([0, 1, 2], dtype='int32')),
             ('offsets', np.array([3], dtype='int32')),
             ('types', np.array([5], dtype='uint8')),
             ('rhoe', np.array([20.0]))],
        ]
        path = write_temp(build_vtu(None, pieces=pieces))
        self.paths.append(path)
        got = vtu.read(path)

        self.assertEqual(len(got['types']), 2, 'a piece was dropped')
        self.assertEqual(len(got['Position']), 6)
        np.testing.assert_array_equal(got['rhoe'], [10.0, 20.0])

    def test_connectivity_is_shifted_into_the_joined_numbering(self):
        """Each piece numbers points from zero, so indices must be rebased.

        Without the shift the second piece's triangle would silently point at
        the first piece's vertices - a coherent-looking mesh made of the wrong
        cells, which is far worse than an error.
        """
        pieces = [
            [('Position', np.zeros((3, 3))),
             ('connectivity', np.array([0, 1, 2], dtype='int32')),
             ('offsets', np.array([3], dtype='int32')),
             ('types', np.array([5], dtype='uint8'))],
            [('Position', np.ones((4, 3))),
             ('connectivity', np.array([0, 1, 2, 3], dtype='int32')),
             ('offsets', np.array([4], dtype='int32')),
             ('types', np.array([5], dtype='uint8'))],
        ]
        path = write_temp(build_vtu(None, pieces=pieces))
        self.paths.append(path)
        got = vtu.read(path)

        # Second piece's indices rebased by the first piece's 3 points.
        np.testing.assert_array_equal(got['connectivity'], [0, 1, 2, 3, 4, 5, 6])
        # offsets are cumulative across the joined connectivity.
        np.testing.assert_array_equal(got['offsets'], [3, 7])
        self.assertLess(int(got['connectivity'].max()), len(got['Position']))

    def test_both_layouts_agree(self):
        """The two layouts carry the same numbers; only the framing differs."""
        wide = self._read(header_type='UInt64', declare=True)
        narrow = self._read(header_type='UInt32', declare=False)
        self.assertEqual(sorted(wide), sorted(narrow))
        for name in wide:
            np.testing.assert_array_equal(wide[name], narrow[name])


if __name__ == '__main__':
    unittest.main(verbosity=2)
