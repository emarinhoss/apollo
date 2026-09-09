"""Minimal reader for the VTK XML files Apollo writes.

Apollo does not write these files itself: it hands a DM to PETSc's
PETSCVIEWERVTK (src/solvers/apsolver.cc:417), so the layout on disk is PETSc's
and changes with the PETSc version and the rank count. Only what the tests need
is implemented: the point coordinates, the cell connectivity, and the data
arrays. Pulling in the VTK bindings to check a number in a test is not worth it.

THE SIZE PREFIX IS NOT A CONSTANT. Each appended data block is preceded by its
length, and the width of that length is whatever the file's `header_type`
attribute says. The VTK default, when the attribute is absent, is UInt32.

    PETSc 3.15   no header_type attribute   4-byte prefix (sizeof(int))
    PETSc 3.19   header_type="UInt64"       8-byte prefix (sizeof(PetscInt64))

Both were read off PETSc's own offset arithmetic in
src/dm/impls/plex/plexvtu.c, where `boffset` advances by the data size plus
`sizeof(int)` at v3.15.5 and `sizeof(PetscInt64)` at v3.19.6. This reader
assumed the 8-byte form unconditionally, which is what Ubuntu 24.04 happens to
ship; on 22.04 - and so in Google Colab - it failed with numpy's "buffer size
must be a multiple of element size", naming neither the file nor the cause.

THERE IS ONE <Piece> PER MPI RANK, AND ALL OF THEM COUNT. PETSc writes a Piece
per rank into the single file, each re-declaring the same array names at its own
offset. Reading DataArray elements into a flat dictionary keeps only the last
rank's copy and silently discards the others - which is where the belief that
"Apollo writes roughly 1/N of the cells on N ranks" came from. The measurements
behind it (7792, 3896, 1961 cells on 1, 2 and 4 ranks) are the size of the last
piece, not of the file: a two-rank run here holds 2152 + 2152 = 4304 cells,
exactly matching the one-rank file. Nothing was ever missing from the output.
"""

import struct
import xml.etree.ElementTree as ET

try:
    import numpy as np
except ImportError:  # pragma: no cover - tests skip themselves without numpy
    np = None

# Base type codes; the byte-order prefix is taken from the file.
_DTYPES = {
    'Float64': 'f8', 'Float32': 'f4',
    'Int64': 'i8', 'Int32': 'i4',
    'UInt64': 'u8', 'UInt32': 'u4', 'UInt8': 'u1',
}

# VTK's header_type, and the struct code for each width. UInt32 is the default
# the VTK format specifies when the attribute is absent, which is the case for
# every file PETSc wrote before 3.16.
_HEADER_TYPES = {'UInt32': 'I', 'UInt64': 'Q'}
_DEFAULT_HEADER_TYPE = 'UInt32'

_APPENDED = b'<AppendedData encoding="raw">'


class _Layout(object):
    """How this particular file frames its appended data blocks."""

    def __init__(self, path, blob, root):
        self.path = path
        self.blob = blob

        # Anything compressed has a different block layout entirely - a header
        # of several integers rather than one - so refuse instead of returning
        # noise.
        compressor = root.get('compressor')
        if compressor:
            raise ValueError('%s: compressed appended data (compressor=%r) is '
                             'not supported by this reader' % (path, compressor))

        self.endian = '>' if root.get('byte_order') == 'BigEndian' else '<'
        self.header_type = root.get('header_type') or _DEFAULT_HEADER_TYPE
        if self.header_type not in _HEADER_TYPES:
            raise ValueError(
                '%s: unsupported header_type=%r (expected one of %s)'
                % (path, self.header_type, ', '.join(sorted(_HEADER_TYPES))))
        self.fmt = self.endian + _HEADER_TYPES[self.header_type]
        self.width = struct.calcsize(self.fmt)

    def block(self, payload, da):
        """One DataArray's values, or None if it is not an appended block."""
        name, dtype, offset = da.get('Name'), da.get('type'), da.get('offset')
        if name is None or offset is None or dtype not in _DTYPES:
            return None

        blob, path = self.blob, self.path
        at = payload + int(offset)
        if at + self.width > len(blob):
            raise ValueError(
                '%s: array %r starts past the end of the file (offset %s, '
                'header_type %s)' % (path, name, offset, self.header_type))
        nbytes = struct.unpack(self.fmt, blob[at:at + self.width])[0]

        # A length that does not fit the file, or is not a whole number of
        # elements, means the prefix was not read where the file put it -
        # almost always the wrong header_type. Say so here, with the numbers,
        # instead of leaving numpy to report a buffer-size mismatch.
        np_dtype = np.dtype(self.endian + _DTYPES[dtype])
        end = at + self.width + nbytes
        if end > len(blob) or nbytes % np_dtype.itemsize:
            raise ValueError(
                '%s: array %r declares %d bytes at offset %s, which is %s '
                '(file is %d bytes, element size %d). The file says '
                'header_type=%s, so the length prefix was read as %d bytes; '
                'if that is wrong every array will be misread.'
                % (path, name, nbytes, offset,
                   'past the end of the file' if end > len(blob)
                   else 'not a whole number of elements',
                   len(blob), np_dtype.itemsize, self.header_type, self.width))

        values = np.frombuffer(blob[at + self.width:end], dtype=np_dtype)
        ncomp = int(da.get('NumberOfComponents', 1))
        return values.reshape(-1, ncomp) if ncomp > 1 else values


def read(path):
    """Return {array name: numpy array} for one .vtu file, all pieces joined."""
    if np is None:
        raise RuntimeError('numpy is required to read VTU files')

    with open(path, 'rb') as fh:
        blob = fh.read()

    marker = blob.find(_APPENDED)
    if marker < 0:
        raise ValueError('%s: no raw appended data section' % path)
    # The payload starts just after the '_' that follows the tag.
    payload = blob.index(b'_', marker) + 1

    # The header is well-formed XML once the unparsed payload is cut away.
    root = ET.fromstring(blob[:marker].decode('utf-8', 'replace') + '</VTKFile>')
    layout = _Layout(path, blob, root)

    # Pieces are concatenated in file order. connectivity indexes points within
    # its own piece and offsets index into its own connectivity, so both are
    # shifted by the running totals as pieces are appended.
    collected = {}
    point_base = 0
    conn_base = 0

    for piece in root.iter('Piece'):
        piece_points = 0
        piece_conn = 0
        for da in piece.iter('DataArray'):
            values = layout.block(payload, da)
            if values is None:
                continue
            name = da.get('Name')
            if name == 'connectivity':
                values = values + point_base
                piece_conn = len(values)
            elif name == 'offsets':
                values = values + conn_base
            elif name == 'Position':
                piece_points = len(values)
            collected.setdefault(name, []).append(values)
        point_base += piece_points
        conn_base += piece_conn

    return {name: parts[0] if len(parts) == 1 else np.concatenate(parts)
            for name, parts in collected.items()}
