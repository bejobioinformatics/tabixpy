import sys

from ._io     import getFilenames, saveVcfGzPy, loadVcfGzPy
from ._tabix  import getPos
from ._logger import logger, getLogLevel

"""
(cd ..; python -c 'import tabixpy; tabixpy.genVcfGzPy("tests/annotated_tomato_150.100000.vcf.gz")')

import tabixpy; _= tabixpy.genVcfGzPy("tests/annotated_tomato_150.100000.vcf.gz")
import tabixpy; _= tabixpy.loadVcfGzPy("tests/annotated_tomato_150.100000.vcf.gz")

import tabixpy; _= tabixpy.genVcfGzPy("tests/annotated_tomato_150.SL2.50ch00-01-02.vcf.gz")
import tabixpy; _= tabixpy.loadVcfGzPy("tests/annotated_tomato_150.SL2.50ch00-01-02.vcf.gz")

import tabixpy; _= tabixpy.genVcfGzPy("tests/annotated_tomato_150.vcf.bgz")
import tabixpy; _= tabixpy.loadVcfGzPy("tests/annotated_tomato_150.vcf.bgz")
"""

def _parseBGZ(filehandle):
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
    ce             = 0

    while block_len >= 0:
        _, first_pos, last_pos, block_len, block_size, chrom_name, num_cols, num_rows = getPos(filehandle, lastReal, 0, None)

        if block_len < 0:
            logger.info(f"block_len {block_len}")

        e  += 1
        ce += 1

        if e == 1 or e % 1_000 == 0:
            logger.info(f"getAllPositions :: reading block {e:12,d} {ce:12,d}")
            logger.info(f"  lastReal {lastReal:12,d} first_pos {first_pos:12,d} last_pos {last_pos:12,d} block_len {block_len:7,d} block_size {block_size:7,d} num_cols {num_cols:5,d} num_rows {num_rows:5,d} chrom_name {chrom_name}")

            # if e%500_000 == 0:
            #     break

        if numCols is None:
            numCols = num_cols

        if chrom_name is None:
            logger.debug(f"chrom_name {chrom_name} is NONE")
            lastReal += block_len
            continue

        if lastChrom != chrom_name:
            if e > 1:
                logger.info(f"getAllPositions :: read    block {e-1:12,d} {ce-1:12,d} {lastChrom}")

            logger.info(f'getAllPositions :: new chrom :: {chrom_name}')
            logger.info(f"  lastReal {lastReal:12,d} first_pos {first_pos:12,d} last_pos {last_pos:12,d} block_len {block_len:7,d} block_size {block_size:7,d} num_cols {num_cols:5,d} num_rows {num_rows:5,d} chrom_name {chrom_name}")

            lastChrom = chrom_name
            ce        = 1
            chroms    .append(chrom_name)
            reals     .append([])
            firsts    .append([])
            lasts     .append([])
            rows      .append([])

        if getLogLevel() == "DEBUG":
            logger.debug(f"getAllPositions {e:4,d} :: lastReal {lastReal:9,d} first_pos {first_pos:9,d} last_pos {last_pos:9,d} block_len {block_len:7,d} block_size {block_size:7,d} chrom_name {chrom_name} num_cols {num_cols:3,d} num_rows {num_rows:3,d}")

        reals [-1].append(lastReal)
        firsts[-1].append(first_pos)
        lasts [-1].append(last_pos)
        rows  [-1].append(num_rows)

        lastReal += block_len

    logger.info(f"getAllPositions :: read    block {e-1:12,d} {ce-1:12,d} {lastChrom}")

    res = {
        "chroms"        : chroms,
        "numCols"       : numCols,
        "chromSizes"    : [len(d) for d in rows],
        "realPositions" : reals,
        "firstPositions": firsts,
        "lastPositions" : lasts,
        "numberRows"    : rows,
    }

    res["chromLength"] = sum(res["chromSizes"])
    
    return res

def readBGZ(infile, save=False):
    # setLogLevel(logging.DEBUG)

    logger.info(f"reading {infile}")

    (ingz, inid, inbj, inbk) = getFilenames(infile)

    inf        = open(ingz, "rb")
    data       = _parseBGZ(inf)

    if getLogLevel() == "DEBUG":
        chroms     = data["chroms"]
        numCols    = data["numCols"]
        reals      = data["realPositions"]
        firsts     = data["firstPositions"]
        lasts      = data["lastPositions"]
        rows       = data["numberRows"]

        logger.debug(f"chroms  {chroms}")
        logger.debug(f"numCols {numCols}")
        logger.debug(f"reals   {reals[0][:10]}-{reals[-1][-10:]}")
        logger.debug(f"firsts  {firsts[0][:10]}-{firsts[-1][-10:]}")
        logger.debug(f"lasts   {lasts[0][:10]}-{lasts[-1][-10:]}")
        logger.debug(f"rows    {rows[0][:10]}-{rows[-1][-10:]}")

    if save:
        saveVcfGzPy(ingz, data)

    inf.close()

    return data

def main(infile):
    genVcfGzPy(infile)
    loadVcfGzPy(infile)

if __name__ == '__main__':
    main(sys.argv[1])
