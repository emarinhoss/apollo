r"""Provides an interface to read FMF files. For documentation on FMF file
format please see http://lanl.arxiv.org/abs/0904.1299v1.
"""

import warnings
warnings.simplefilter('ignore', UserWarning)

class FMFSection:
    r"""Represents a generic section in the FMF file. [*data definition] and
    [*data] are treated differently.
    """

    def __init__(self, name):
        r"""__init__(str:name) -> FMFSection object.

        Create a new generic FMF section object with given name."""

        # name of section
        self.name = name
        # data stored as key-value pairs
        self.data = {}

    def getKeys(self):
        return self.data.keys()

    def getValue(self, name):
        return self.data[name]
