import gzip
import json
import struct

from ._consts import COMPRESS, GZIP_MAGIC, TABIXPY_FORMAT_NAME, TABIXPY_FORMAT_VER
from ._logger     import logger, getLogLevel

def genStructValueGetter(fhd, returnBytes=False):
    def getValues(fmt):
        fmt_s = struct.calcsize(fmt)
        pack  = fhd.read(fmt_s)
        res   = struct.unpack(fmt, pack)
        if returnBytes:
            return res, pack
        else:
            return res
    return getValues

def getFilenames(infile, old_index_ext=".tbi", new_index_ext=".tbj"):
    if infile.endswith(old_index_ext):
        ingz     = infile[:-4]
        inid     = infile
        inbj     = ingz + new_index_ext
    else:
        ingz     = infile
        inid     = infile + old_index_ext
        inbj     = ingz + new_index_ext

    return ingz, inid, inbj

def load(ingz, format_name=TABIXPY_FORMAT_NAME, format_ver=TABIXPY_FORMAT_VER):
    _, _, inbj = getFilenames(ingz)

    compressed = None
    with open(inbj, "rb") as fhd:
        firstChars = fhd.read(2)
        if firstChars[0] == 123: # { 123
            compressed = False
        else:
            assert firstChars == GZIP_MAGIC, firstChars
            compressed = True

    is_compressed = "compressed json" if compressed else "json"
    logger.info(f"loading {inbj} as {is_compressed}")

    data = None
    if compressed:
        with gzip.open(inbj, "rb") as fhd:
            data = json.load(fhd)
    else:
        with open(inbj, "rt") as fhd:
            data = json.load(fhd)

    assert "__format_name__" in data
    assert "__format_ver__"  in data

    assert data["__format_name__"] == format_name
    assert data["__format_ver__" ] == format_ver

    return data

def save(data, ingz, compress=COMPRESS, ext=".tbj", format_name=TABIXPY_FORMAT_NAME, format_ver=TABIXPY_FORMAT_VER):
    data["__format_name__"] = format_name
    data["__format_ver__" ] = format_ver

    outfileJ = ingz + ext

    logger.info(f"saving  {outfileJ}")

    opener   = open
    if compress:
        logger.debug("compressing")
        opener = gzip.open

    with opener(outfileJ, "wt") as fhd:
        json.dump(data, fhd, indent=1)
