import os
import sys
import json
import struct
import hashlib

from ._io     import genStructValueGetter, getFilenames
from ._gzip   import gzip, GZIP_MAGIC
from ._tabix  import getPos
from ._logger import logger, getLogLevel
from ._consts import (
    COMPRESS,
    VCFBGZ_FORMAT_VER,
    VCFBGZ_FORMAT_NAME,
    VCFBGZ_EXTENSION,
    VCFBGZ_EOF
)

"""
python
(cd ..; python -c 'import tabixpy; tabixpy.genAllBlocks("tests/annotated_tomato_150.100000.vcf.gz")')
import tabixpy; tabixpy.genAllBlocks("tests/annotated_tomato_150.100000.vcf.gz")
import tabixpy; _= tabixpy.loadAllBlocks("tests/annotated_tomato_150.100000.vcf.gz")

import tabixpy; tabixpy.genAllBlocks("tests/annotated_tomato_150.SL2.50ch00-01-02.vcf.gz")
import tabixpy; _= tabixpy.loadAllBlocks("tests/annotated_tomato_150.SL2.50ch00-01-02.vcf.gz")
"""

def getAllBlocks(filehandle):
    block_len      = 0
    lastReal       = 0
    chroms         = []
    lastChrom      = None
    numCols        = None
    reals          = []
    firsts         = []
    lasts          = []
    rows           = []
    e              = 0

    while block_len >= 0:
        _, first_pos, last_pos, block_len, block_size, chrom_name, num_cols, num_rows = getPos(filehandle, lastReal, 0, None)
        
        if block_len < 0:
            break
        
        if numCols is None:
            numCols = num_cols

        if lastChrom != chrom_name:
            lastChrom = chrom_name
            chroms    .append(chrom_name)
            reals     .append([])
            firsts    .append([])
            lasts     .append([])
            rows      .append([])

        if getLogLevel() == "DEBUG":
            logger.debug(f"getAllPositions {e:4,d} :: lastReal {lastReal:9,d} first_pos {first_pos:9,d} last_pos {last_pos:9,d} block_len {block_len:9,d} block_size {block_size:9,d} chrom_name {chrom_name} num_cols {num_cols:3,d} num_rows {num_rows:3,d}")
        
        reals [-1].append(lastReal)
        firsts[-1].append(first_pos)
        lasts [-1].append(last_pos)
        rows  [-1].append(num_rows)

        lastReal += block_len
        e        += 1

    res = {
        "chroms": chroms,
        "numCols": numCols,
        "chromSizes": [len(d) for d in rows],
        "realPositions": reals,
        "firstPositions": firsts,
        "lastPositions": lasts,
        "numberRows": rows,
    }

    res["chromLength"] = sum(res["chromSizes"])
    
    return res

def saveAllBlocks(filename, data, compress=COMPRESS, ext=VCFBGZ_EXTENSION, format_name=VCFBGZ_FORMAT_NAME, format_ver=VCFBGZ_FORMAT_VER):
    logger.info(f" saving {filename}{ext}")

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
        ["q"                   , len(format_name)     ],
        [f"{len(format_name)}s", format_name.encode() ],
        ["q"                   , format_ver           ],
        ["q"                   , len(headerJ)         ],
        [f"{len(headerJ)}s"    , headerJ.encode()     ]
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

    with opener(filename + ext, 'wb') as fhd:
        fhd.write(header_dat)

        for lstK in ["realPositions", "firstPositions", "lastPositions", "numberRows"]:
            logger.info(f" writing {lstK} - {chromLength}")
            lst  = data[lstK]
            
            for chrom_data in lst:
                # logger.info(f"chrom_data {chrom_data[:10]} {chrom_data[-10:]}")
                st = chrom_data[0]
                chrom_data = [st] + [v - chrom_data[c] for c,v in enumerate(chrom_data[1:])]
                # logger.info(f"chrom_data {chrom_data[:10]} {chrom_data[-10:]}")
                lstD = struct.pack(f"<{len(chrom_data)}q"     , *chrom_data)
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

def loadAllBlocks(filename, ext=VCFBGZ_EXTENSION, format_name=VCFBGZ_FORMAT_NAME, format_ver=VCFBGZ_FORMAT_VER):
    indexFile = filename + ext
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
                raise ValueError(f"not a valid uncompressed file. invalid magic header: {fmt}. expected {GZIP_MAGIC} OR {format_name}")
            
            if fmt == VCFBGZ_FORMAT_NAME:
                compressed = False
            else:
                raise ValueError(f"not a valid uncompressed file. invalid magic header: {fmt}. expected {GZIP_MAGIC} OR {format_name}")

        if compressed is None:
            raise ValueError(f"not a valid uncompressed file. invalid magic header: {fmt}. expected {GZIP_MAGIC} OR {format_name}")

    opener   = open
    if compressed:
        logger.info(" decompressing")
        opener = gzip.open

    with opener(filename + ext, 'rb') as fhd:
        getter      = genStructValueGetter(fhd, returnBytes=True)
        
        ((fmt_len, ), d) = getter("<q")
        m.update(d)
        logger.debug(f" fmt_len    {fmt_len}")

        assert fmt_len == len(format_name), f"fmt_len {fmt_len} == len(format_name) {len(format_name)}"

        ((fmt_nam, ), d) = getter(f"<{fmt_len}s")
        m.update(d)
        fmt_nam     = fmt_nam.decode()
        logger.debug(f" fmt_nam    {fmt_nam}")

        assert fmt_nam == format_name, f"fmt_nam {fmt_nam} == format_name {format_name}"

        ((fmt_ver, ), d) = getter("<q")
        m.update(d)
        logger.debug(f" fmt_ver    {fmt_ver}")
        assert fmt_ver == format_ver, f"fmt_ver {fmt_ver} == format_ver {format_ver}"

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
                logger.info(f"  reading {chromSize} values")
                (chrom_data, d) = getter(f"<{chromSize}q")
                m.update(d)

                chrom_data = list(chrom_data)
                # logger.info(f"chrom_data {chrom_data[:10]} {chrom_data[-10:]}")
                for c in range(1,len(chrom_data)):
                    chrom_data[c] = chrom_data[c] + chrom_data[c-1]
                # logger.info(f"chrom_data {chrom_data[:10]} {chrom_data[-10:]}")

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

def genAllBlocks(infile):
    # setLogLevel(logging.DEBUG)

    logger.info(f"reading {infile}")

    ingz, inid, inbj = getFilenames(infile)

    inf        = open(ingz, "rb")
    get_values = genStructValueGetter(inf)
    data       = getAllBlocks(inf)

    chroms     = data["chroms"]
    numCols    = data["numCols"]
    reals      = data["realPositions"]
    firsts     = data["firstPositions"]
    lasts      = data["lastPositions"]
    rows       = data["numberRows"]

    if getLogLevel() == "DEBUG":
        logger.debug(f"chroms  {chroms}")
        logger.debug(f"numCols {numCols}")
        logger.debug(f"reals   {reals[0][:10]}-{reals[-1][-10:]}")
        logger.debug(f"firsts  {firsts[0][:10]}-{firsts[-1][-10:]}")
        logger.debug(f"lasts   {lasts[0][:10]}-{lasts[-1][-10:]}")
        logger.debug(f"rows    {rows[0][:10]}-{rows[-1][-10:]}")

    saveAllBlocks(ingz, data)
    # return data

def main(infile):
    pass

if __name__ == '__main__':
    main(sys.argv[1])