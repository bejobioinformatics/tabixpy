import os
import sys
import json

from ._gzip       import openGzipStream
from ._io         import getFilenames, loadTabixPy, saveTabixPy
from ._logger     import logger, setLogLevel
from ._tabix      import readTabix

class Tabix:
    def __init__(self, ingz, logLevel=None):
        self._infile  = ingz
        self._numCols = None

        if logLevel is not None:
            setLogLevel(logLevel)

        self._ingz, self._inid, self._inbj, self._inbk = getFilenames(self._infile)

        if os.path.exists(self._inbj):
            logger.info("reading TBJ file")
            self.load()

        else:
            logger.info("reading Tabix file")
            self._data   = readTabix(self._ingz)

    def load(self):
        self._data = loadTabixPy(self._ingz)

    def save(self, overwrite=True, compress=True):
        if os.path.exists(self._inbj):
            if not overwrite:
                return
        saveTabixPy(self._data, self._ingz, compress=compress)
    
    @property
    def data(self):
        return self._data

    @property
    def chromosomes(self):
        return self._data["names"]

    @property
    def numCols(self):
        return self._numCols

    def getChromosomeIter(self, chrom, begin=None, end=None, asLine=False):
        assert  chrom in self.chromosomes
        assert  begin is None or begin >= 0
        assert  end   is None or end   >= 0
        assert (begin is None or end is None) or begin <= end

        numCols = self.numCols

        idx     = self.chromosomes.index(chrom)

        ref     = self._data["refs"   ][idx]
        intvs   = ref["intvs"]

        first_block = ref["first_block"]
        last_block  = ref["last_block" ]

        if "block_len" not in first_block:
            raise NotImplementedError("block search not implemented")

        assert ref["ref_name"] == chrom
        assert end is None or (end <= last_block["last_pos"]), f"(end {end:12,d} <= last_block['last_pos'] {last_block['last_pos']:12,d}) - last_block {last_block}"

        intvsBegin = None
        if begin is not None:
            for rpos, r in enumerate(intvs):
                if r["first_pos"] == begin:
                    r["chunk"] = rpos
                    intvsBegin = r
                    break
                elif r["first_pos"] > begin:
                    if rpos == 0: #first block
                        r["chunk_n"] = rpos
                        intvsBegin = r
                        break
                    else:
                        r = intvs[rpos - 1] #went too far. go back one block
                        r["chunk_n"] = rpos - 1
                        intvsBegin = r
                        break
            if intvsBegin is None: #not in any block. use last block
                r = intvs[-1]
                r["chunk_n"] = -1
                intvsBegin = r # last block
        else:
            r = intvs[0]
            r["chunk_n"] = 0
            intvsBegin = r

        # logger.debug(f"begin      {begin}")
        # logger.debug(f"end        {end}"  )
        # logger.debug(f"intvs[ 0]  {intvs[ 0]}")
        # logger.debug(f"intvs[-2]  {intvs[-2]}")
        # logger.debug(f"intvs[-1]  {intvs[-1]}")
        logger.debug(f"intvsBegin {intvsBegin}")

        with openGzipStream(self._ingz, intvsBegin["real"], 0, asLine=asLine, chrom=chrom, begin=begin, end=end) as fhd:
            for line in fhd:
                yield line


def main(infile):
    tabix  = Tabix(infile)
    tabix.save(overwrite=True)

if __name__ == "__main__":
    infile = sys.argv[1]
    main(infile)
