import sys

from ._io     import getFilenames, saveVcfGzPy, loadVcfGzPy
from ._tabix  import getPos
from ._logger import logger, getLogLevel

"""
python
(cd ..; python -c 'import tabixpy; tabixpy.genVcfGzPy("tests/annotated_tomato_150.100000.vcf.gz")')
import tabixpy; tabixpy.genVcfGzPy("tests/annotated_tomato_150.100000.vcf.gz")
import tabixpy; _= tabixpy.loadVcfGzPy("tests/annotated_tomato_150.100000.vcf.gz")

import tabixpy; tabixpy.genVcfGzPy("tests/annotated_tomato_150.SL2.50ch00-01-02.vcf.gz")
import tabixpy; _= tabixpy.loadVcfGzPy("tests/annotated_tomato_150.SL2.50ch00-01-02.vcf.gz")
"""

def getVcfGzPy(filehandle):
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

def genVcfGzPy(infile):
    # setLogLevel(logging.DEBUG)

    logger.info(f"reading {infile}")

    ingz, _, _ = getFilenames(infile)

    inf        = open(ingz, "rb")
    data       = getVcfGzPy(inf)

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

    saveVcfGzPy(ingz, data)
    # return data

def main(infile):
    genVcfGzPy(infile)
    loadVcfGzPy(infile)

if __name__ == '__main__':
    main(sys.argv[1])