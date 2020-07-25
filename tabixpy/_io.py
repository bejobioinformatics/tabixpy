import gzip
import json
import struct
import hashlib

from ._logger     import logger, getLogLevel
from ._gzip       import GZIP_MAGIC

from ._consts import (
    COMPRESS,
)

from ._consts import (
    TABIX_EXTENSION,
)

from ._consts import (
    TABIXPY_FORMAT_NAME,
    TABIXPY_FORMAT_VER,
    TABIXPY_EXTENSION
)

from ._consts import (
    VCFBGZ_FORMAT_VER,
    VCFBGZ_FORMAT_NAME,
    VCFBGZ_EXTENSION,
    VCFBGZ_EOF
)


BIN_SIZES = [
    [      0, 2** 8, 'B' ], # char               1
    [ -2** 7, 2** 7, 'b' ], # char               1

    [      0, 2**16, 'H' ], # short unsigned     2
    [ -2**15, 2**15, 'h' ], # short              2

    [      0, 2**32, 'L' ], # long unsigned      4
    [ -2**31, 2**31, 'l' ], # long               4

    [      0, 2**64, 'Q' ], # long long unsigned 8
    [ -2**63, 2**63, 'q' ], # long long          8
]

def getByteSize(vals):
    for (min_val, max_val, code) in BIN_SIZES:
        if all([x >= min_val and x < max_val for x in vals]):
            return code
    
    raise ValueError(f"not able to encode values {vals}")

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

def getFilenames(infile):
    if infile[-3:] == '.gz' or infile[-4:] == '.bgz':
        ingz     = infile
    else:
        ingz     = infile[:-4]

    inid     = infile + TABIX_EXTENSION
    inbj     = ingz   + TABIXPY_EXTENSION
    inbk     = ingz   + VCFBGZ_EXTENSION

    assert not ingz.endswith(".")

    return ingz, inid, inbj, inbk

def saveTabixPy(ingz, data, compress=COMPRESS):
    data["__format_name__"] = TABIXPY_FORMAT_NAME
    data["__format_ver__" ] = TABIXPY_FORMAT_VER

    outfileJ = ingz + TABIXPY_EXTENSION

    logger.info(f"saving  {outfileJ}")

    opener   = open
    if compress:
        logger.debug("compressing")
        opener = gzip.open

    with opener(outfileJ, "wt") as fhd:
        json.dump(data, fhd, indent=1)

def loadTabixPy(ingz):
    (ingz, inid, inbj, inbk) = getFilenames(ingz)

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

    assert data["__format_name__"] == TABIXPY_FORMAT_NAME
    assert data["__format_ver__" ] == TABIXPY_FORMAT_VER

    return data

def saveVcfGzPy(filename, data, compress=COMPRESS):
    outfile = filename + VCFBGZ_EXTENSION

    logger.info(f" saving {outfile}")

    flatten     = lambda lst: (item for sublist in lst for item in sublist)
    chromLength = data["chromLength"]
    header      = {
        "chroms"     : data["chroms"],
        "numCols"    : data["numCols"],
        "chromSizes" : data["chromSizes"],
        "chromLength": data["chromLength"]
    }

    headerJ     = json.dumps(header)

    header_fmts = [
        ["q"                          , len(VCFBGZ_FORMAT_NAME)     ],
        [f"{len(VCFBGZ_FORMAT_NAME)}s", VCFBGZ_FORMAT_NAME.encode() ],
        ["q"                          , VCFBGZ_FORMAT_VER           ],
        ["q"                          , len(headerJ)                ],
        [f"{len(headerJ)}s"           , headerJ.encode()            ]
    ]
    
    
    m = hashlib.sha256()

    header_fmt  = "<" + "".join([h[0] for h in header_fmts])
    logger.debug(f"header_fmt '{header_fmt}'")
    header_val  = [h[1] for h in header_fmts]
    header_dat  = struct.pack(header_fmt, *header_val )
    logger.debug(header_dat)
    logger.debug(header_val)
    m.update(header_dat)

    if getLogLevel() == "DEBUG":
        header_rev  = struct.unpack(header_fmt, header_dat)
        header_rev  = list(header_rev)
        logger.debug(header_rev)
        assert header_rev == header_val

    opener   = open
    if compress:
        logger.info(" compressing")
        opener = gzip.open

    with opener(outfile, 'wb') as fhd:
        fhd.write(header_dat)

        for lstK in ["realPositions", "firstPositions", "lastPositions", "numberRows"]:
            logger.info(f" writing {lstK:16s} - {chromLength:18,d}")
            lst  = data[lstK]
            
            for chrom_data in lst:
                # logger.info(f"chrom_data {chrom_data[:10]} {chrom_data[-10:]}")
                cdsum      = sum(chrom_data)
                st         = chrom_data[0]
                chrom_data = [st] + [v - chrom_data[c] for c,v in enumerate(chrom_data[1:])]
                # logger.info(f"chrom_data {chrom_data[:10]} {chrom_data[-10:]}")
                
                fmt = getByteSize(chrom_data)
                fms = f"<qc{len(chrom_data)}{fmt}"
                logger.info(f"   fmt {fmt} min {min(chrom_data):15,d} max {max(chrom_data):15,d} len {len(chrom_data):18,d} cdsum {cdsum:21,d} fmts {fms}")

                lstD       = struct.pack(fms, cdsum, fmt.encode(), *chrom_data)
                # lstD       = struct.pack(f"<{len(chrom_data)}q"     , *chrom_data)

                fhd.write(lstD)
                m.update(lstD)
                # sys.exit(0)

            # lstD = struct.pack(f"<{chromLength}q"     , *flatten(lst))
            # fhd.write(lstD)
            # m.update(lstD)

        digestHex  = m.hexdigest()
        digestLen  = len(digestHex)
        digestSize = struct.pack(f"<q", digestLen)
        m.update(digestSize)
        fhd.write(digestSize)
        digestHex  = m.hexdigest()

        digest = struct.pack(f"<{digestLen}s", digestHex.encode())
        fhd.write(digest)
        
        logger.info(digestHex)

        fhd.write(VCFBGZ_EOF)

    return 

def loadVcfGzPy(filename):
    indexFile = filename + VCFBGZ_EXTENSION
    logger.info(f" loading {indexFile}")

    m = hashlib.sha256()

    compressed = None
    with open(indexFile, "rb") as fhd:
        firstChars = fhd.read( 8 + len(VCFBGZ_FORMAT_NAME) )
        compressed = None

        if firstChars[:2] == GZIP_MAGIC:
            compressed = True
        else:
            fmt = firstChars[8:]
            
            try:
                fmt = fmt.decode()
            except:
                raise ValueError(f"not a valid uncompressed file. invalid magic header: {fmt}. expected {GZIP_MAGIC} OR {VCFBGZ_FORMAT_NAME}")
            
            if fmt == VCFBGZ_FORMAT_NAME:
                compressed = False
            else:
                raise ValueError(f"not a valid uncompressed file. invalid magic header: {fmt}. expected {GZIP_MAGIC} OR {VCFBGZ_FORMAT_NAME}")

        if compressed is None:
            raise ValueError(f"not a valid uncompressed file. invalid magic header: {fmt}. expected {GZIP_MAGIC} OR {VCFBGZ_FORMAT_NAME}")

    opener   = open
    if compressed:
        logger.info(" decompressing")
        opener = gzip.open

    with opener(indexFile, 'rb') as fhd:
        getter      = genStructValueGetter(fhd, returnBytes=True)
        
        ((fmt_len, ), d) = getter("<q")
        m.update(d)
        logger.debug(f" fmt_len    {fmt_len}")

        assert fmt_len == len(VCFBGZ_FORMAT_NAME), f"fmt_len {fmt_len} == len(VCFBGZ_FORMAT_NAME) {len(VCFBGZ_FORMAT_NAME)}"

        ((fmt_nam, ), d) = getter(f"<{fmt_len}s")
        m.update(d)
        fmt_nam     = fmt_nam.decode()
        logger.debug(f" fmt_nam    {fmt_nam}")

        assert fmt_nam == VCFBGZ_FORMAT_NAME, f"fmt_nam {fmt_nam} == VCFBGZ_FORMAT_NAME {VCFBGZ_FORMAT_NAME}"

        ((fmt_ver, ), d) = getter("<q")
        m.update(d)
        logger.debug(f" fmt_ver    {fmt_ver}")
        assert fmt_ver == VCFBGZ_FORMAT_VER, f"fmt_ver {fmt_ver} == VCFBGZ_FORMAT_VER {VCFBGZ_FORMAT_VER}"

        ((lenHeaderJ, ), d) = getter("<q")
        m.update(d)
        logger.debug(f" lenHeaderJ {lenHeaderJ}")

        ((headerJ, ), d) = getter(f"<{lenHeaderJ}s")
        m.update(d)
        headerJ     = headerJ.decode()
        header      = json.loads(headerJ)
        logger.debug(f" header     {header}")

        chromLength = header["chromLength"]

        for lstK in ["realPositions", "firstPositions", "lastPositions", "numberRows"]:
            logger.info(f" reading {lstK}")
            header[lstK] = []
            for chromSize in header["chromSizes"]:
                logger.info(f"  {chromSize:12,d} values")

                ((cdsum,), d) = getter(f"<q")
                m.update(d)

                ((fmt,), d) = getter(f"<c")
                m.update(d)
                # logger.info(f"cdsum {cdsum} fmt {fmt}")

                (chrom_data, d) = getter(f"<{chromSize}{fmt.decode()}")
                m.update(d)

                chrom_data = list(chrom_data)
                # logger.info(f"chrom_data {chrom_data[:10]} {chrom_data[-10:]}")
                for c in range(1,len(chrom_data)):
                    chrom_data[c] = chrom_data[c] + chrom_data[c-1]
                # logger.info(f"chrom_data {chrom_data[:10]} {chrom_data[-10:]}")

                if cdsum == 0: #pypy
                    assert sum(chrom_data) == -1, f"sum(chrom_data) {sum(chrom_data)} == cdsum {cdsum}"

                else:
                    assert sum(chrom_data) == cdsum, f"sum(chrom_data) {sum(chrom_data)} == cdsum {cdsum}"

                header[lstK].append(chrom_data)

        ((digestLen, ), d) = getter("<q")
        m.update(d)
        logger.debug(f"digestLen  {digestLen}")

        ((digestHex, ), _)  = getter(f"<{digestLen}s")
        digestHex = digestHex.decode()
        logger.info(f"digestHex  {digestHex}")
        assert digestHex == m.hexdigest()

        eof = fhd.read(len(VCFBGZ_EOF))
        
        assert eof == VCFBGZ_EOF

        assert len(fhd.read()) == 0

    return header
