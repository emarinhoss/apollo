"""Minimal reader for the VTK XML files Apollo writes.

Apollo does not write these files itself: it hands a DM to PETSc's
PETSCVIEWERVTK (src/solvers/apsolver.cc:417), so the layout on disk is PETSc's
and **changes with the PETSc version**. Only what the tests need is
implemented: the point coordinates, the cell connectivity, and the cell-data
arrays. Pulling in the VTK bindings to check a number in a test is not worth it.

THE SIZE PREFIX IS NOT A CONSTANT. Each appended data block is preceded by its
length, and the width of that length is whatever the file's `header_type`
attribute says. The VTK default, when the attribute is absent, is UInt32.

    PETSc 3.15   no header_type attribute   4-byte prefix (sizeof(int))
    PETSc 3.19   header_type="UInt64"       8-byte prefix (sizeof(PetscInt64))

Both were read off PETSc's own offset arithmetic in
src/dm/impls/plex/plexvtu.c, where `boffset` advances by the data size plus
`sizeof(int)` at v3.15.5 and `sizeof(PetscInt64)` at v3.19.6.

This reader assumed the 8-byte form unconditionally, which is what Ubuntu 24.04
happens to ship. On Ubuntu 22.04 - and so in Google Colab - it read four bytes
of real data as part of the length and failed with numpy's

    ValueError: buffer size must be a multiple of element size

naming neither the file nor the cause. It now reads the attribute, and when a
block still does not make sense it says which file, which array and what it
tried, rather than letting numpy guess.
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


def read(path):
    """Return {array name: numpy array} for one .vtu file."""
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

    # Anything compressed has a different block layout entirely - a header of
    # several integers rather than one - so refuse instead of returning noise.
    compressor = root.get('compressor')
    if compressor:
        raise ValueError('%s: compressed appended data (compressor=%r) is not '
                         'supported by this reader' % (path, compressor))

    endian = '>' if root.get('byte_order') == 'BigEndian' else '<'

    header_type = root.get('header_type') or _DEFAULT_HEADER_TYPE
    if header_type not in _HEADER_TYPES:
        raise ValueError('%s: unsupported header_type=%r (expected one of %s)'
                         % (path, header_type, ', '.join(sorted(_HEADER_TYPES))))
    header_fmt = endian + _HEADER_TYPES[header_type]
    header_width = struct.calcsize(header_fmt)

    arrays = {}
    for da in root.iter('DataArray'):
        name, dtype, offset = da.get('Name'), da.get('type'), da.get('offset')
        if name is None or offset is None or dtype not in _DTYPES:
            continue
        at = payload + int(offset)
        if at + header_width > len(blob):
            raise ValueError(
                '%s: array %r starts past the end of the file (offset %s, '
                'header_type %s)' % (path, name, offset, header_type))
        nbytes = struct.unpack(header_fmt, blob[at:at + header_width])[0]

        # A length that does not fit the file, or is not a whole number of
        # elements, means the prefix was not read where the file put it -
        # almost always the wrong header_type. Say so here, with the numbers,
        # instead of leaving numpy to report a buffer-size mismatch.
        np_dtype = np.dtype(endian + _DTYPES[dtype])
        end = at + header_width + nbytes
        if end > len(blob) or nbytes % np_dtype.itemsize:
            raise ValueError(
                '%s: array %r declares %d bytes at offset %s, which is %s '
                '(file is %d bytes, element size %d). The file says '
                'header_type=%s, so the length prefix was read as %d bytes; '
                'if that is wrong every array will be misread.'
                % (path, name, nbytes, offset,
                   'past the end of the file' if end > len(blob)
                   else 'not a whole number of elements',
                   len(blob), np_dtype.itemsize, header_type, header_width))

        values = np.frombuffer(blob[at + header_width:end], dtype=np_dtype)
        ncomp = int(da.get('NumberOfComponents', 1))
        arrays[name] = values.reshape(-1, ncomp) if ncomp > 1 else values
    return arrays
