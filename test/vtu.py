"""Minimal reader for the VTK XML files Apollo writes.

Apollo emits UnstructuredGrid files with raw appended data and a UInt64 byte
count in front of each block. That is a documented VTK layout, but no
general-purpose Python package in this repository's dependency list reads it,
and pulling in the VTK bindings just to check a number in a test is not worth
it. Only what the tests need is implemented: the point coordinates, the cell
connectivity, and the cell-data arrays.
"""

import struct
import xml.etree.ElementTree as ET

try:
    import numpy as np
except ImportError:  # pragma: no cover - tests skip themselves without numpy
    np = None

_DTYPES = {
    'Float64': '<f8', 'Float32': '<f4',
    'Int64': '<i8', 'Int32': '<i4',
    'UInt64': '<u8', 'UInt32': '<u4', 'UInt8': '<u1',
}

_APPENDED = b'<AppendedData encoding="raw">'


def read(path):
    """Return {array name: numpy array} for one .vtu file."""
    if np is None:
        raise RuntimeError('numpy is required to read VTU files')

    with open(path, 'rb') as fh:
        blob = fh.read()

    marker = blob.find(_APPENDED)
    if marker < 0:
        raise ValueError(f'{path}: no raw appended data section')
    # The payload starts just after the '_' that follows the tag.
    payload = blob.index(b'_', marker) + 1

    # The header is well-formed XML once the unparsed payload is cut away.
    root = ET.fromstring(blob[:marker].decode('utf-8', 'replace') + '</VTKFile>')

    arrays = {}
    for da in root.iter('DataArray'):
        name, dtype, offset = da.get('Name'), da.get('type'), da.get('offset')
        if name is None or offset is None or dtype not in _DTYPES:
            continue
        at = payload + int(offset)
        nbytes = struct.unpack('<Q', blob[at:at + 8])[0]
        values = np.frombuffer(blob[at + 8:at + 8 + nbytes], dtype=_DTYPES[dtype])
        ncomp = int(da.get('NumberOfComponents', 1))
        arrays[name] = values.reshape(-1, ncomp) if ncomp > 1 else values
    return arrays
