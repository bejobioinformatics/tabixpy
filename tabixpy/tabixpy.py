import os
import sys
import json
import bisect

from ._gzip       import openGzipStream
from ._io         import getFilenames
from ._io         import loadTabixPy, saveTabixPy
from ._io         import loadVcfGzPy, saveVcfGzPy
from ._logger     import logger, setLogLevel, getLogLevel
from ._tabix      import readTabix
from ._vcfbgzpy   import readBGZ

from enum import Enum, auto
class Formats(Enum):
    TBJ       = auto()
    TBK       = auto()
    DEFAULT   = TBJ

class TabixDefaults():
    OVERWRITE = True
    COMPRESS  = True
    CINE      = True

class Tabix(TabixDefaults):
    def __init__(self, ingz, indexType=Formats.DEFAULT, logLevel=None):
        self._infile    = ingz
        self._indexType = indexType
        self._numCols   = None
        self._data      = None
        self._type      = None

        if logLevel is not None:
            setLogLevel(logLevel)

        self._inbgz, self._intbi, self._intbj, self._intbk = getFilenames(self._infile)

    def create(self, overwrite=TabixDefaults.OVERWRITE, compress=TabixDefaults.COMPRESS):
        if   self._indexType == Formats.TBJ:
            self.createTBJ(overwrite=overwrite, compress=compress)
        elif self._indexType == Formats.TBK:
            self.createTBK(overwrite=overwrite, compress=compress)
        else:
            raise ValueError(f"NO SUCH INDEX TYPE {self._indexType}. Valid values are TBJ and TBK")

    def createTBJ(self, overwrite=TabixDefaults.OVERWRITE, compress=TabixDefaults.COMPRESS):
        logger.info(f"creating TBJ")
        if self._data is None or self._type != Formats.TBJ:
            self.loadTBI()
        self.saveTBJ(overwrite=overwrite, compress=compress)

    def createTBK(self, overwrite=TabixDefaults.OVERWRITE, compress=TabixDefaults.COMPRESS):
        logger.info(f"creating TBK")
        if self._data is None or self._type != Formats.TBK:
            self.loadBGZ()
        self.saveTBK(overwrite=overwrite, compress=compress)

    def save(self, overwrite=TabixDefaults.OVERWRITE, compress=TabixDefaults.COMPRESS):
        if   self._type == Formats.TBJ:
            self.saveTBJ(overwrite=overwrite, compress=compress)

        elif self._type == Formats.TBK:
            self.saveTBK(overwrite=overwrite, compress=compress)

        else:
            raise ValueError(f"NO SUCH INDEX TYPE {self._type}. Valid values are DEFAULT (TBK), TBJ and TBK")

    def saveTBJ(self, overwrite=TabixDefaults.OVERWRITE, compress=TabixDefaults.COMPRESS):
        logger.info(f"saving TBJ")

        if self._type != Formats.TBJ:
            raise ValueError("Saving TBJ when data is in different format")

        if os.path.exists(self._intbj):
            if not overwrite:
                return

        saveTabixPy(self._inbgz, self._data, compress=compress)

    def saveTBK(self, overwrite=TabixDefaults.OVERWRITE, compress=TabixDefaults.COMPRESS):
        logger.info(f"saving TBK")

        if self._type != Formats.TBK:
            raise ValueError("Saving TBK when data is in different format")

        if os.path.exists(self._intbk):
            if not overwrite:
                return

        saveVcfGzPy(self._inbgz, self._data)

    def load(self, create_if_not_exists=TabixDefaults.CINE):
        if self._indexType == Formats.TBK:
            if   os.path.exists(self._intbk) or create_if_not_exists:
                logger.info(f"reading TBK index {self._intbk}")
                self.loadTBK(create_if_not_exists=create_if_not_exists)
            else:
                raise IOError(f"no such index {self._intbk} OR FILE {self._inbgz}")

        elif self._indexType == Formats.TBJ:
            if   os.path.exists(self._intbj) or create_if_not_exists:
                logger.info(f"reading TBJ index {self._intbj}")
                self.loadTBJ(create_if_not_exists=create_if_not_exists)
            else:
                raise IOError(f"no such TBJ index {self._intbj} or TBI index {self._intbi}")

        else:
            raise ValueError(f"NO SUCH INDEX TYPE {self._type}. Valid values are DEFAULT (TBK), TBJ and TBK")

    def loadBGZ(self, create_if_not_exists=TabixDefaults.CINE):
        if not os.path.exists(self._inbgz):
            raise IOError(f"BGZ file {self._inbgz} does not exists")

        logger.info(f"loading BGZ")

        self._data = readBGZ(self._inbgz)
        self._type = Formats.TBK

    def loadTBI(self, create_if_not_exists=TabixDefaults.CINE):
        if not os.path.exists(self._intbi):
            raise IOError(f"TBI file {self._intbi} does not exists")

        logger.info(f"loading TBI")

        self._data  = readTabix(self._inbgz)
        self._type = Formats.TBJ

    def loadTBJ(self, create_if_not_exists=TabixDefaults.CINE):
        if not os.path.exists(self._intbj):
            if not create_if_not_exists:
                raise IOError(f"TBJ file {self._intbj} does not exists")
            else:
                self.loadTBI()
                self.saveTBJ()

        logger.info(f"loading TBJ")

        self._data = loadTabixPy(self._inbgz)
        self._type = Formats.TBJ

    def loadTBK(self, create_if_not_exists=TabixDefaults.CINE):
        if not os.path.exists(self._intbk):
            if not create_if_not_exists:
                raise IOError(f"TBK file {self._intbk} does not exists")
            else:
                self.loadBGZ()
                self.saveTBK()

        logger.info(f"loading TBK")

        self._data = loadVcfGzPy(self._inbgz)
        self._type = Formats.TBK
    
    @property
    def bgz(self):
        return self._inbgz

    @property
    def indexFile(self):
        if self._indexType == Formats.TBJ:
            return self._intbj

        if self._indexType == Formats.TBK:
            return self._intbk

    @property
    def sourceFile(self):
        if self._indexType == Formats.TBJ:
            return self._intbi

        if self._indexType == Formats.TBK:
            return self._inbgz

    @property
    def data(self):
        return self._data

    @property
    def chromosomes(self):
        if self._indexType == Formats.TBJ:
            return self._data["names"]

        if self._indexType == Formats.TBK:
            return self._data["chroms"]

    @property
    def numCols(self):
        return self._data.get("numCols", None)

    def getChromosomeIter(self, chrom, begin=None, end=None, asLine=False):
        assert  chrom in self.chromosomes
        assert  begin is None or begin >= 0
        assert  end   is None or end   >= 0
        assert (begin is None or end is None) or begin <= end

        if self._indexType == Formats.TBJ:
            return self.getChromosomeIterTBJ(chrom, begin=begin, end=end, asLine=asLine)

        if self._indexType == Formats.TBK:
            return self.getChromosomeIterTBK(chrom, begin=begin, end=end, asLine=asLine)

    def getChromosomeIterTBJ(self, chrom, begin=None, end=None, asLine=False):
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

        with openGzipStream(self._inbgz, intvsBegin["real"], 0, asLine=asLine, chrom=chrom, begin=begin, end=end) as fhd:
            for line in fhd:
                yield line

    def getChromosomeIterTBK(self, chrom, begin=None, end=None, asLine=False):
        numCols    = self.numCols

        idx        = self.chromosomes.index(chrom)

        reals      = self._data["realPositions"][idx]
        firsts     = self._data["firstPositions"][idx]
        lasts      = self._data["lastPositions"][idx]
        rows       = self._data["numberRows"][idx]

        t2l        = lambda x: ",".join([f"{y:12,d}" for y in x])

        if getLogLevel() == "DEBUG":
            logger.debug(f"getChromosomeIterTBK :: chrom {chrom} begin {begin if begin is not None else -1:12,d} end {end if end is not None else -1:12,d} asLine {asLine}")
            logger.debug(f"getChromosomeIterTBK ::   numCols {numCols}")
            logger.debug( "getChromosomeIterTBK ::   reals   {} | {}".format(t2l(reals[:4] ), t2l(reals[-4:]) ))
            logger.debug( "getChromosomeIterTBK ::   firsts  {} | {}".format(t2l(firsts[:4]), t2l(firsts[-4:])))
            logger.debug( "getChromosomeIterTBK ::   lasts   {} | {}".format(t2l(lasts[:4] ), t2l(lasts[-4:]) ))
            logger.debug( "getChromosomeIterTBK ::   rows    {} | {}".format(t2l(rows[:4]  ), t2l(rows[-4:])  ))

        if begin is not None and begin > lasts[-1]:
            logger.debug(f"getChromosomeIterTBK :: EMPTY :: begin {begin if begin is not None else -1:12,d} <= lasts[-1] {lasts[-1]:12,d}")
            return []

        if end   is not None and end   > lasts[-1]:
            logger.debug(f"getChromosomeIterTBK :: EMPTY :: end   {end   if end   is not None else -1:12,d} <= lasts[-1] {lasts[-1]:12,d}")
            return []

        real = 0
        if begin is not None:
            pos  = bisect.bisect_left(firsts, begin)

            if pos >= len(firsts):
                pos = len(firsts) - 1

            posb   = pos - 1 if pos - 1 >= 0          else pos
            posa   = pos + 1 if pos + 1 < len(firsts) else pos

            firstb = firsts[posb]
            first  = firsts[pos]
            firsta = firsts[posa]

            real   = reals[pos ]
            realb  = reals[posb]

            logger.debug(f"POSITION begin {begin if begin is not None else -1:12,d} - {posb:12,d} ({firstb:12,d}) < {pos:12,d} ({first:12,d}) < {posa:12,d} ({firsta:12,d}) - {len(firsts):12,d}")

            if first > begin:
                real = realb

            if firstb > begin:
                real = -1

        if real == -1:
            if idx == 0:
                raise ValueError(f"begin {begin} smallest then begining {firsts[0]}")
            else:
                logger.debug("reverting to previous chromosome")
                real = self._data["realPositions"][idx - 1][-1]

        logger.debug(f"POSITION :: real {real:12,d}")

        with openGzipStream(self._inbgz, real, 0, asLine=asLine, chrom=chrom, begin=begin, end=end) as fhd:
            for line in fhd:
                yield line
