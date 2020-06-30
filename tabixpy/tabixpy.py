import os
import sys
import gzip
import struct
import json
import logging

__format_ver__    = 5
__format_name__   = "TBJ"

COMPRESS          = True
BRICK_SIZE        = 216
BLOCK_SIZE        = 2**16
FILE_BYTES_MASK   = 0xFFFFFFFFFFFFF0000
BLOCK_BYTES_MASK  = 0x0000000000000FFFF
MAX_BIN           = (((1<<18)-1)//7)

# Custom formatter
#https://stackoverflow.com/questions/14844970/modifying-logging-message-format-based-on-message-logging-level-in-python3
class CustomConsoleFormatter(logging.Formatter):
    _fmts = {
        # logging.NOSET   : "NOSET   : %(msg)s",
        logging.DEBUG   : "DEBUG   : %(asctime)8s - %(name)-7s/%(module)-7s/%(funcName)-12s/%(lineno)3d - %(msg)s",
        logging.INFO    : "INFO    : %(msg)s",
        logging.WARNING : "WARNING : %(msg)s",
        logging.ERROR   : "ERROR   : %(name)-7s/%(module)-7ss/%(funcName)-12s/%(lineno)3d - %(msg)s",
        logging.CRITICAL: "CRITICAL: %(name)-7s/%(module)-7ss/%(funcName)-12s/%(lineno)3d - %(msg)s",
    }
    # format  = "%(asctime)8s - %(name)-7s/%(funcName)-12s/%(lineno)3d - %(levelname)-7s - %(message)s",
    # datefmt = "%H:%M:%S",
    # datefmt='%Y-%m-%d %H:%M:%S',

    def __init__(self, fmt="%(levelno)d: %(msg)s", datefmt="%H:%M:%S"):
        super().__init__(fmt=fmt, datefmt=datefmt, style='%')

    def format(self, record):
        # Save the original format configured by the user
        # when the logger formatter was instantiated
        format_orig = self._style._fmt

        # Replace the original format with one customized by logging level
        self._style._fmt = self._fmts.get(record.levelno, format_orig)

        # Call the original formatter class to do the grunt work
        result = logging.Formatter.format(self, record)

        # Restore the original format configured by the user
        self._style._fmt = format_orig

        return result


# Set up a logger
logger          = logging.getLogger('tabixpy')

my_formatter    = CustomConsoleFormatter(datefmt="%H:%M:%S",)

console_handler = logging.StreamHandler()
console_handler.setFormatter(my_formatter)

logger.addHandler(console_handler)
logger.setLevel(logging.INFO)



class Tabix:
    def __init__(self, ingz, logLevel=None):
        self._infile  = ingz
        self._numCols = None

        if logLevel is not None:
            setLogLevel(logLevel)

        self._ingz, self._inid, self._inbj = get_filenames(self._infile)

        if os.path.exists(self._inbj):
            logger.info("reading TBJ file")
            self._data  = load(self._ingz)

        else:
            logger.info("reading Tabix file")
            self._data   = read_tabix(self._ingz)

    def save(self, overwrite=True, compress=True):
        if os.path.exists(self._inbj):
            if not overwrite:
                return
        save(self._data, self._ingz, compress=compress)
    
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

class openGzipStream():
    def __init__(self, infile, realPos, bytesPos, asLine=False, chrom=None, begin=None, end=None):
        logger.debug(f"infile {infile} realPos {realPos} bytesPos {bytesPos} asLine {asLine} chrom {chrom} begin {begin} end {end}")

        self._infile    = infile
        self._realPos   = realPos
        self._bytesPos  = bytesPos
        self._asLine    = asLine
        self._chrom     = chrom
        self._begin     = begin
        self._end       = end
        self._fhdf      = None
        self._fhdz      = None

    def __enter__(self):
        self._fhdf = open(self._infile, 'rb')
        
        currPos = self._realPos
        self._fhdf.seek(self._realPos)
        
        self._fhdz = gzip.open(self._fhdf, 'rb')

        if self._bytesPos is not None and self._bytesPos > 0:
            raise NotImplementedError("BytesPos Not implemented")
            # block    = block[self._bytesPos:]
            self._fhdz.seek(self._bytesPos)
            block    = self._fhdz.read(BLOCK_SIZE-self._bytesPos).decode()
        else:
            block    = self._fhdz.read(BLOCK_SIZE).decode()


        foundChrom = False
        numCols    = None
        lastLine   = None
        prevLine   = None
        prevCols   = None
        rn         = -1
        bn         = -1
        while len(block) > 0:
            bn       += 1
            lines     = block.split("\n")
            block     = self._fhdz.read(BLOCK_SIZE).decode()

            if lastLine is not None:
                lines[0] = lastLine + lines[0]
                lastLine = None

            logger.debug(f"bn {bn} len(lines) {len(lines)} :: filter")
            lines   = [line for line in lines if len(line) > 0 and line[0] != "#"] # filter out empty and comment lines 
            columns = [line.split("\t") for line in lines]

            if len(columns[-1]) != len(columns[1]): #incomplete last line
                logger.debug(f"bn {bn} len(lines) {len(lines)} :: pop")
                lastLine = lines.pop()
                columns.pop()

            if len(columns[0]) != len(columns[1]): #incomplete first line
                logger.debug(f"bn {bn} len(lines) {len(lines)} :: shift")
                lines   = lines[1:]
                columns = columns[1:]

            if numCols is None:
                numCols = len(columns[0])

            if self._chrom is not None:
                logger.debug(f"bn {bn} len(lines) {len(lines)} :: chrom {self._chrom}")
                lines   = [line for ln, line in enumerate(lines) if columns[ln][0] == self._chrom]
                columns = [col  for col      in columns          if col[0]         == self._chrom]

            if len(lines) == 0:
                logger.debug(f"bn {bn} len(lines) == 0 :: return")
                return

            for ln, line in enumerate(lines):
                cols = columns[ln]
                logger.debug(f"bn {bn} len(lines) == {len(lines)} :: chrom {self._chrom} begin {self._begin} end {self._end} :: ln {ln}")

                if self._begin is not None or self._end is not None:
                    pos = int(cols[1])
                    
                    if self._begin is not None and pos < self._begin:
                        logger.debug(f" pos {pos} < self._begin {self._begin}")
                        continue

                    if self._end is not None and pos >= self._end:
                        logger.debug(f"pos {pos} >= self._end {self._end}")
                        return

                    if self._asLine:
                        yield line
                    else:
                        yield cols

                else:
                    if self._asLine:
                        yield line
                    else:
                        yield cols


    def __exit__(self, type, value, traceback):
        self._fhdz.close()
        self._fhdf.close()

# samtools format specs:
# https://samtools.github.io/hts-specs/SAMv1.pdf
# https://github.com/xbrianh/bgzip/blob/master/bgzip/__init__.py
bgzip_eof = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")

class EOF:
    __slot__ = []
    pass

def read_gzip_header(fp):
    #https://github.com/python/cpython/blob/3.8/Lib/gzip.py
    
    #http://www.htslib.org/doc/bgzip.html
    # the gzip header includes an extra sub-field with identifier 'BC' 
    # and the length of the compressed block, including all headers.
    #
    # http://samtools.github.io/hts-specs/SAMv1.pdf
    # 1. The F.EXTRA bit in the header is set to indicate that extra fields
    #    are present.
    #
    # 2. The extra field used by BGZF uses the two subfield ID values 66 and 67
    #    (ascii ‘BC’)
    #
    # 3. The length of the BGZF extra field payload (field LEN in the gzip
    #    specification) is 2 (two bytes of payload).
    #
    # 4. The payload of the BGZF extra field is a 16-bit unsigned integer in
    #    little endian format. This integer gives the size of the containing
    #    BGZF block minus one.

    
    FTEXT, FHCRC, FEXTRA, FNAME, FCOMMENT = 1, 2, 4, 8, 16

    tell  = fp.tell()

    magic = fp.read(2)
    if magic == b'':
        return False

    if magic != b'\037\213':
        raise IOError('Not a gzipped file (%r)' % magic)

    (method, flag, last_mtime) = struct.unpack("<BBIxx", fp.read(8))

    if method != 8:
        raise IOError('Unknown compression method')

    if flag & FEXTRA:
        # Read & discard the extra field, if present
        extra_len, = struct.unpack("<H", fp.read(2))
        extra_data = fp.read(extra_len)
        # logger.debug(f" extra length {extra_len}")
        # logger.debug(f" extra data   {extra_data}")
        if extra_len >= 4:
            if b"BC" in extra_data:
                bci            = extra_data.index(b'BC')
                subfield_len_d = extra_data[bci+2:bci+2+2]
                subfield_len,  = struct.unpack("<H", subfield_len_d)
                # logger.debug(f" subfieldlen   {subfield_len} {subfield_len_d}")
                
                block_len_d    = extra_data[bci+2+2:bci+2+2+subfield_len]
                block_len,     = struct.unpack("<H", block_len_d)
                # logger.debug(f" block len   {block_len + 1:12,d}")

                fp.seek(tell)

                return block_len + 1
    
    fp.seek(tell)

    return None

def setLogLevel(level):
    logger.setLevel(level)

def getLogLevel():
    return logging.getLevelName(logger.getEffectiveLevel())

def gen_value_getter(fhd):
    def get_values(fmt):
        fmt_s = struct.calcsize(fmt)
        res   = struct.unpack(fmt, fhd.read(fmt_s))
        return res
    return get_values

def get_filenames(infile):
    if infile.endswith(".tbi"):
        ingz     = infile[:-4]
        inid     = infile
        inbj     = ingz + ".tbj"
    else:
        ingz     = infile
        inid     = infile + ".tbi"
        inbj     = ingz + ".tbj"

    return ingz, inid, inbj

def load(ingz):
    _, _, inbj = get_filenames(ingz)

    compressed = None
    with open(inbj, "rb") as fhd:
        firstChars = fhd.read(2)
        if firstChars[0] == 123: # { 123
            compressed = False
        else:
            assert firstChars == b"\x1F\x8B", firstChars
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

    assert data["__format_name__"] == __format_name__
    assert data["__format_ver__" ] == __format_ver__

    return data

def save(data, ingz, compress=COMPRESS):
    data["__format_name__"] = __format_name__
    data["__format_ver__" ] = __format_ver__

    outfileJ = ingz + ".tbj"

    logger.info(f"saving  {outfileJ}")

    opener   = open
    if compress:
        logger.debug("compressing")
        opener = gzip.open

    with opener(outfileJ, "wt") as fhd:
        json.dump(data, fhd, indent=1)

def reg2bin(begPos, endPos):
    #- Given a zero-based, half-closed and half-open region [beg, end), 
    # the bin number is calculated with the following C function:

    endPos -= 1
    if (begPos>>14 == endPos>>14): return ((1<<15)-1)//7 + (begPos>>14)
    if (begPos>>17 == endPos>>17): return ((1<<12)-1)//7 + (begPos>>17)
    if (begPos>>20 == endPos>>20): return ((1<< 9)-1)//7 + (begPos>>20)
    if (begPos>>23 == endPos>>23): return ((1<< 6)-1)//7 + (begPos>>23)
    if (begPos>>26 == endPos>>26): return ((1<< 3)-1)//7 + (begPos>>26)
    return 0

def reg2bins(rbeg, rend):
    # The list of bins that may overlap a region [beg, end) can be obtained 
    # with the following C function:
    #define MAX_BIN (((1<<18)-1)/7)
    
    res = [None] * MAX_BIN
    
    i       = 0
    k       = 0
    rend   -= 1
    res[i]  = 0
    i      += 1

    ranges = [
        [    1 + (rbeg>>26),      1 + (rend>>26)],
        [    9 + (rbeg>>23),      9 + (rend>>23)],
        [   73 + (rbeg>>20),     73 + (rend>>20)],
        [  585 + (rbeg>>17),    585 + (rend>>17)],
        [ 4681 + (rbeg>>14),   4681 + (rend>>14)]
    ]

    for b,e in ranges:
        for k in range(b ,e+1):
            res[i] = k
            i+=1

    return i; # #elements in list[]

def getBlock(filehandle, real_pos, block_len=None):
    filehandle.seek(real_pos, 0)

    if block_len is None:
        block_len = read_gzip_header(filehandle)

    if block_len is None:
        raise IOError("invalid block length")

    if block_len == 0:
        raise IOError("block length is zero")

    filehandle.seek(real_pos, 0)

    # logger.debug(f" block len   {block_len:12,d} {block_len//1024:12,d}kb {block_len//1024//1024:12,d}mb")

    blockdata = filehandle.read(block_len)
    block     = gzip.decompress(blockdata)
    blocktext = block.decode()

    assert len(blockdata) == block_len
    assert len(block) == len(blocktext)

    # logging.error(f"len(block) {len(blockdata)} == header_len {header_len} + len(bgzip_eof) {len(bgzip_eof)}")
    if len(blockdata) == len(bgzip_eof):
        if blockdata == bgzip_eof:
            return EOF(), -1

    # logger.debug(f" header_len  {header_len:12,d}")
    # logger.debug(f" blockdata   {len(blockdata):12,d}")
    # logger.debug(f" block       {len(block):12,d}")
    # logger.debug(f" blocktext   {len(blocktext):12,d}")

    assert len(blockdata) == block_len, f"block data has wrong size :: real_pos {real_pos} block_len {block_len}"
    assert len(blockdata) > 0         , f"block data has size zero :: real_pos {real_pos} block_len {block_len}"
    assert len(block) > 0             , f"block has size zero :: real_pos {real_pos} block_len {block_len}"
    assert len(blocktext) > 0         , f"block text has size zero :: real_pos {real_pos} block_len {block_len}"

    # rows = block.decode().split("\n")

    # logging.info(len(rows))
    # logging.info(rows[0])
    # logging.info(rows[-2])
    # logging.info(rows[-1])

    return blocktext, block_len

def parseBlock(block, bytes_pos, chrom):
    bin_pos   = -1
    first_pos = -1
    last_pos  = -1

    assert len(block) > 0, f"empty block '{block}'"
    assert bytes_pos  < len(block)

    rows      = [row for row in block.split("\n") if len(row) > 0 and row[0] != "#"]
    first_row = None
    last_row  = None
    rows      = [row for row in rows if row.split("\t")[0] == chrom]
    num_rows  = len(rows)

    if num_rows == 0:
        # logger.debug(f"num_rows {num_rows} == 0")
        return -1, -1, -1, -1, -1, -1

    elif num_rows == 1:
        first_row = rows[0]
        last_row  = rows[0]
        # logger.debug(f"num_rows {num_rows} == 1 :: first_row {first_row[:2]} last_row {last_row[:2]}")
    
    else:
        first_cols  = rows[ 0].split("\t")
        second_cols = rows[ 1].split("\t")

        if len(first_cols) != len(second_cols) or first_cols[0] != second_cols[0]:
            first_cols  = second_cols
            num_rows   -= 1
        
        assert first_cols[0] == chrom
        
        for i in range(len(rows)):
            # logger.debug(f"i {i} (i*-1)-1 {(i*-1)-1:3d} len(rows) {len(rows):3d}")
            last_cols   = rows[(i*-1)-1].split("\t")
            if len(last_cols) >= 2 and last_cols[0] == first_cols[0] and last_cols[0] == chrom and int(last_cols[1]) > int(first_cols[1]):
                # logger.debug(f"i {i} (i*-1)-1 {(i*-1)-1:3d} len(rows) {len(rows):3d} num_rows {num_rows:3d}\n\tlen(last_cols) {len(last_cols):3d} == len(first_cols) {len(first_cols):3d} or last_cols[0] {last_cols[0]} == first_cols[0] {first_cols[0]} - last_cols[1] {last_cols[1]} == first_cols[1] {first_cols[1]}\n")
                break
            else:
                # logger.debug(f"i {i} (i*-1)-1 {(i*-1)-1:3d} len(rows) {len(rows):3d} num_rows {num_rows:3d}\n\tlen(last_cols) {len(last_cols):3d} != len(first_cols) {len(first_cols):3d} or last_cols[0] {last_cols[0]} != first_cols[0] {first_cols[0]} - last_cols[1] {last_cols[1]} == first_cols[1] {first_cols[1]}")
                num_rows   -= 1

    assert len(first_cols) >= len(last_cols), f"{len(first_cols)} == {len(last_cols)}\n{len(first_cols)} {first_cols}\n{len(last_cols)} {last_cols}"

    first_chrom = first_cols[0]
    first_pos   = first_cols[1]
    first_pos   = int(first_pos)

    last_chrom  = last_cols[0]
    last_pos    = last_cols[1]
    last_pos    = int(last_pos)

    assert first_chrom == last_chrom, f"first_chrom {first_chrom} == last_chrom {last_chrom}"

    bin_reg    = block[bytes_pos:].split("\n")[0]
    bin_cols   = bin_reg.split("\t")
    if len(bin_cols) == len(first_cols):
        bin_chrom  = bin_cols[0]
        bin_pos    = bin_cols[1]
        
        assert len(bin_cols) == len(first_cols), f"len(bin_cols) {len(bin_cols)} == len(first_cols) {len(first_cols)} :: bytes_pos {bytes_pos} chrom {chrom}\nblock[         :100] {block[:100]}\nblock[bytes_pos:   ] {block[bytes_pos:bytes_pos+100]}"

        try:
            bin_pos    = int(bin_pos)
        except ValueError as e:
            logging.error(e)
            logging.error(bin_reg)
            raise

        assert len(bin_cols  ) >= len(last_cols), f"{len(bin_cols  )} == {len(last_cols)}\n{len(bin_cols  )} {bin_cols}\n{len(last_cols)} {last_cols}"
            
    return bin_pos, first_pos, last_pos, first_chrom, len(first_cols), num_rows

def getPos(filehandle, real_pos, bytes_pos, chrom):
    # logger.debug(f"getPos :: real_pos {real_pos} bytes_pos {bytes_pos} chrom {chrom}")

    bin_pos    = -1
    first_pos  = -1
    last_pos   = -1
    block_len  = -1
    block_size = -1
    chrom_name = None
    num_cols   = -1
    num_rows   = -1

    block, block_len = getBlock(filehandle, real_pos)

    if isinstance(block, EOF):
        return bin_pos, first_pos, last_pos, block_len, block_size, chrom_name, num_cols, num_rows 

    assert len(block) > 0, f"got empty block {real_pos}"

    block_size       = len(block)

    bin_pos, first_pos, last_pos, chrom_name, num_cols, num_rows = parseBlock(block, bytes_pos, chrom)

    return bin_pos, first_pos, last_pos, block_len, block_size, chrom_name, num_cols, num_rows

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
        _, first_pos, last_pos, block_len, block_size, chrom_name, num_cols, num_rows = getPos(filehandle, lastReal, 0)
        
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

        logger.debug(f"getAllPositions {e:4,d} :: lastReal {lastReal:9,d} first_pos {first_pos:9,d} last_pos {last_pos:9,d} block_len {block_len:9,d} block_size {block_size:9,d} chrom_name {chrom_name} num_cols {num_cols:3,d} num_rows {num_rows:3,d}")
        
        reals [-1].append(lastReal)
        firsts[-1].append(first_pos)
        lasts [-1].append(last_pos)
        rows  [-1].append(num_rows)

        lastReal += block_len
        e        += 1

    return {
        "chroms": chroms,
        "numCols": numCols,
        "realPositions": reals,
        "firstPositions": firsts,
        "lastPositions": lasts,
        "numberRows": rows
    }

def gen_all_blocks(infile):
    logger.info(f"reading {infile}")

    ingz, inid, inbj = get_filenames(infile)

    inf        = open(ingz, "rb")
    get_values = gen_value_getter(fhd)
    data       = getAllBlocks(inf)

def read_tabix(infile):
    logger.info(f"reading {infile}")

    ingz, inid, inbj = get_filenames(infile)

    assert os.path.exists(inid), inid

    fhd        = gzip.open(inid, "rb")
    inf        = open(ingz, "rb")
    get_values = gen_value_getter(fhd)
    data       = {}

    (magic,)   = get_values('<4s')

    assert magic == b'TBI\x01', f'invalid magic: {magic}'

    logger.debug(f"magic                    Magic string                                    char[4]  TBI\1")
    logger.debug(f"-------------------------------------------------------------------------------------------------")

    (
        n_ref,   #                # sequences                                     int32_t
        f_format,#                Format (0: generic; 1: SAM; 2: VCF)             int32_t
        col_seq, #                Column for the sequence name                    int32_t
        col_beg, #                Column for the start of a region                int32_t
        col_end, #                Column for the end of a region                  int32_t
        meta,    #                Leading character for comment lines             int32_t
        skip,    #                # lines to skip at the beginning                int32_t
        l_nm     #                Length of concatenated sequence names           int32_t
    ) = get_values('<8i')

    logger.debug(f"n_ref                    # sequences                                     int32_t  {n_ref:15,d}"   )
    logger.debug(f"format                   Format (0: generic; 1: SAM; 2: VCF)             int32_t  {f_format:15,d}")
    logger.debug(f"col_seq                  Column for the sequence name                    int32_t  {col_seq:15,d}" )
    logger.debug(f"col_beg                  Column for the start of a region                int32_t  {col_beg:15,d}" )
    logger.debug(f"col_end                  Column for the end of a region                  int32_t  {col_end:15,d}" )
    logger.debug(f"meta                     Leading character for comment lines             int32_t  {meta:15,d}"    )
    logger.debug(f"skip                     # lines to skip at the beginning                int32_t  {skip:15,d}"    )
    logger.debug(f"l_nm                     Length of concatenated sequence names           int32_t  {l_nm:15,d}"    )

    assert n_ref   >      0, n_ref
    assert n_ref   <  2**31, n_ref
    assert f_format in [1,2], f_format
    assert col_seq >=     0, col_seq
    assert col_seq <     10, col_seq
    assert col_beg >=     0, col_beg
    assert col_beg <     10, col_beg
    assert col_end >=     0, col_end
    assert col_end <     10, col_end
    assert meta    >      0, meta
    assert meta    <    256, meta
    assert skip    >=     0, skip
    assert skip    <  2**20, skip
    assert l_nm    >      0, l_nm
    assert l_nm    <  2**20, l_nm

    data["n_ref"          ] = n_ref
    data["format"         ] = f_format
    data["col_seq"        ] = col_seq
    data["col_beg"        ] = col_beg
    data["col_end"        ] = col_end
    data["meta"           ] = chr(meta)
    data["skip"           ] = skip
    data["l_nm"           ] = l_nm

    names_s    = fhd.read(l_nm)
    assert len(names_s) == l_nm, len(names_s)

    names = [n.decode("utf-8") for n in names_s.split(b'\x00') if n]
    assert len(names) == n_ref, len(names)

    data["names"  ] = names

    logger.debug(f"names                    Concatenated names, each zero terminated        char[l_nm] [{len(names)}]")
    logger.debug("  - " + "\n  - ".join(names))

    data["refs"   ] = [None] * n_ref

    logger.debug(f"======================== List of indices (n=n_ref [{n_ref:15,d}]) ============================")
    for ref_n in range(n_ref):
        ref                 = {}
        ref["ref_n"   ]     = ref_n
        ref["ref_name"]     = names[ref_n]
        data["refs"][ref_n] = ref

        (n_bin,)            = get_values('<i')
        
        logger.debug(f"     ref_n {ref_n+1:15,d}/{n_ref:15,d} ({names[ref_n]})")
        logger.debug(f"       n_bin             # distinct bins (for the binning index)         int32_t  {n_bin:15,d}")

        assert n_bin >     0, n_bin
        assert n_bin < 2**31, n_bin

        ref["n_bin"     ] = n_bin
        ref["bins"      ] = [None] * n_bin
        ref["bins_begin"] = None
        ref["bins_end"  ] = None

        logger.debug(f"======================== List of distinct bins (n=n_bin [{n_bin:15,d}]) ======================")
        last_chk_real_beg = None
        last_chk_real_end = None
        
        position_memoize = {}
        bin_min_pos      = 2**64
        bin_max_pos      = 0

        for bin_n in range(n_bin):
            (bin_v,n_chunk) = get_values('<Ii')

            logger.debug(f"         bin_n {bin_n+1}/{n_bin}")
            logger.debug(f"           bin           Distinct bin number                             uint32_t {bin_v:15,d}")
            logger.debug(f"           n_chunk       # chunks                                        int32_t  {n_chunk:15,d}")
            
            assert bin_v   >     0, bin_v
            assert bin_v   < 2**31, bin_v
            assert n_chunk >     0, n_chunk
            assert n_chunk < 2**31, n_chunk

            ref["bins"][bin_n] = {
                "bin_n"  : bin_n,
                "bin"    : bin_v,
                "n_chunk": n_chunk,
            }

            chunks = get_values('<' + ('Q'*(n_chunk*2)))

            assert all([i >=     0 for i in chunks]), chunks
            assert all([i <  2**63 for i in chunks])

            # In the compressed file, each uncompressed byte in the text data file is
            # assigned a unique 64-bit virtual file offset where the higher 48 bits keep the
            # real file offset of the start of the gzip block the byte falls in, and the
            # lower 16 bits store the offset of the byte inside the gzip block. Given a
            # virtual file offset, one can directly seek to the start of the gzip block using
            # the higher 48 bits, decompress the block and retrieve the byte pointed by the
            # lower 16 bits of the virtual offset.

            # real_file_offsets   = [c & FILE_BYTES_MASK   for c in chunks]
            real_file_offsets   = [c >> 16              for c in chunks]
            block_bytes_offsets = [c & BLOCK_BYTES_MASK for c in chunks]

            assert all([i >=     0     for i in real_file_offsets])
            assert all([i <  2**48     for i in real_file_offsets])

            assert all([i >=          0 for i in block_bytes_offsets])
            assert all([i <  BLOCK_SIZE for i in block_bytes_offsets])

            chunks_data = {
                "chunk_begin" : [None] * n_chunk,
                "chunk_end"   : [None] * n_chunk,
            }

            for chunk_n, chunk_i in enumerate(range(0,len(chunks),2)):
                chk_beg       = chunks[chunk_i+0]
                chk_end       = chunks[chunk_i+1]

                chk_real_beg  = real_file_offsets[chunk_i+0]
                chk_real_end  = real_file_offsets[chunk_i+1]

                chk_bytes_beg = block_bytes_offsets[chunk_i+0]
                chk_bytes_end = block_bytes_offsets[chunk_i+1]

                # logger.debug(f"begin")
                # logger.debug(f"  block             {chk_beg:064b} {chk_beg:15,d}")
                # logger.debug(f"  mask              {BLOCK_BYTES_MASK:064b}")
                # logger.debug(f"  real offset       {chk_real_beg:064b} {chk_real_beg:15,d}")
                # logger.debug(f"  block byte offset {chk_bytes_beg:064b} {chk_bytes_beg:15,d}")
                # logger.debug(f"end")
                # logger.debug(f"  block             {chk_end:064b} {chk_end:15,d}")
                # logger.debug(f"  mask              {BLOCK_BYTES_MASK:064b}")
                # logger.debug(f"  real offset       {chk_real_end:064b} {chk_real_end:15,d}")
                # logger.debug(f"  block byte offset {chk_bytes_end:064b} {chk_bytes_end:15,d}")
                # logger.debug("")

                assert chk_beg       <  chk_end
                assert chk_real_beg  <= chk_real_end, f"chk_real_beg {chk_real_beg:12,d}  < chk_real_end {chk_real_end:12,d}"
                # assert chk_bytes_beg < chk_bytes_end

                # if last_chk_real_beg is not None:
                #     assert last_chk_real_beg <=  chk_real_beg, f"last_chk_real_beg {last_chk_real_beg:12,d} <  chk_real_beg {chk_real_beg:12,d}"
                #     # if last_chk_real_end is not None:
                #     #     assert last_chk_real_end == chk_real_beg, f"last_chk_real_end {last_chk_real_end:12,d} == chk_real_beg {chk_real_beg:12,d}"
                    
                # if last_chk_real_end is not None:
                #     assert last_chk_real_end <=  chk_real_end
                #     # assert last_chk_real_end == chk_real_beg

                # last_chk_real_beg = chk_real_beg
                # last_chk_real_end = chk_real_end

                if chk_beg in position_memoize:
                    chunk_nfo_beg     = position_memoize[chk_beg]
                    chk_block_len_beg = chunk_nfo_beg["block_len"]
                    chk_bin_pos_beg   = chunk_nfo_beg["bin_pos"]
                    chk_first_pos_beg = chunk_nfo_beg["first_pos"]
                    chk_last_pos_beg  = chunk_nfo_beg["last_pos"] 

                else:
                    chk_bin_pos_beg, chk_first_pos_beg, chk_last_pos_beg, chk_block_len_beg, _, _, _, _ = getPos(inf, chk_real_beg, chk_bytes_beg, ref['ref_name'])
                    
                    chunk_nfo_beg = {
                        "bin_n": bin_n,
                        "chunk_n": chunk_n,
                        "real": chk_real_beg,
                        "bytes": chk_bytes_beg,
                        "block_len": chk_block_len_beg,
                        "bin_pos": chk_bin_pos_beg,
                        "first_pos": chk_first_pos_beg,
                        "last_pos": chk_last_pos_beg,
                    }
                    
                    position_memoize[chk_beg] = chunk_nfo_beg

                if chk_end in position_memoize:
                    chunk_nfo_end     = position_memoize[chk_end]
                    chk_block_len_end = chunk_nfo_end["block_len"]
                    chk_bin_pos_end   = chunk_nfo_end["bin_pos"]
                    chk_first_pos_end = chunk_nfo_end["first_pos"]
                    chk_last_pos_end  = chunk_nfo_end["last_pos"] 

                else:
                    chk_bin_pos_end, chk_first_pos_end, chk_last_pos_end, chk_block_len_end, _, _, _, _ = getPos(inf, chk_real_end, chk_bytes_end, ref['ref_name'])
                
                    chunk_nfo_end = {
                        "bin_n": bin_n,
                        "chunk_n": chunk_n,
                        "real": chk_real_end,
                        "bytes": chk_bytes_end,
                        "block_len": chk_block_len_end,
                        "bin_pos": chk_bin_pos_end,
                        "first_pos": chk_first_pos_end,
                        "last_pos": chk_last_pos_end
                    }

                    position_memoize[chk_end] = chunk_nfo_end

                chunks_data["chunk_begin" ][chunk_n] = chunk_nfo_beg
                
                chunks_data["chunk_end"   ][chunk_n] = chunk_nfo_end

                if chunk_nfo_beg["first_pos"] < bin_min_pos:
                    ref["bins_begin"] = chunks_data["chunk_begin" ][chunk_n]
                    bin_min_pos = chunk_nfo_beg["first_pos"]

                if chunk_nfo_end["last_pos"] > bin_max_pos:
                    ref["bins_end"  ] = chunks_data["chunk_end"   ][chunk_n]
                    bin_max_pos = chunk_nfo_end["last_pos"]

                if getLogLevel() == "DEBUG":
                    logger.debug(f"======================== List of chunks (n=n_chunk[{chunk_n+1:15,d}]) ============================")
                    logger.debug(f"             chunk_n {chunk_n+1:15,d}/{n_chunk:15,d}")
                    logger.debug(f"               cnk_beg   Virtual file offset of the start of the chunk   uint64_t {chk_real_beg:15,d} {chk_bytes_beg:15,d}")
                    logger.debug(f"               cnk_end   Virtual file offset of the end of the chunk     uint64_t {chk_real_end:15,d} {chk_bytes_end:15,d}")
                    logger.debug(f"               cnk_beg_block_len                                         uint64_t {chk_block_len_beg:15,d}")
                    logger.debug(f"               cnk_end_block_len                                         uint64_t {chk_block_len_end:15,d}")
                    logger.debug(f"               cnk_beg_1st_pos                                           uint64_t {chk_first_pos_beg:15,d}")
                    logger.debug(f"               cnk_end_1st_pos                                           uint64_t {chk_first_pos_end:15,d}")
                    logger.debug(f"               cnk_beg_bin_pos                                           uint64_t {chk_bin_pos_beg:15,d}")
                    logger.debug(f"               cnk_end_bin_pos                                           uint64_t {chk_bin_pos_end:15,d}")
                    logger.debug(f"               cnk_beg_lst_pos                                           uint64_t {chk_last_pos_beg:15,d}")
                    logger.debug(f"               cnk_end_lst_pos                                           uint64_t {chk_last_pos_end:15,d}")

            ref["bins"      ][bin_n]["chunks"      ] = chunks_data
            ref["bins"      ][bin_n]["chunks_begin"] = chunks_data["chunk_begin"][ 0]
            ref["bins"      ][bin_n]["chunks_end"  ] = chunks_data["chunk_end"  ][-1]

        # logger.debug(f'bins_begin {ref["bins_begin"]}')
        # logger.debug(f'bins_end   {ref["bins_end"  ]}')

        last_bin = ref["bins_end"  ]
        last_bin_pos, last_first_pos, last_last_pos, last_block_len = [last_bin[k] for k in ['bin_pos', 'first_pos', 'last_pos', 'block_len']]
        logger.debug(f"LAST 0 last_bin_pos {last_bin_pos:12,d} last_first_pos {last_first_pos:12,d} last_last_pos {last_last_pos:12,d} last_block_len {last_block_len:12,d}")

        ref["first_block"] = ref["bins_begin"]
        ref["last_block" ] = ref["bins_end"  ]

        lastReal     = last_bin["real"] + last_block_len
        chrom_name   = ref["ref_name"]
        extra_blocks = []
        e = 0
        while last_block_len > 0 and chrom_name == ref["ref_name"]:
            last_bin_pos, last_first_pos, last_last_pos, last_block_len,          _, chrom_name, _, _ = getPos(inf, lastReal, 0, ref['ref_name'])
            if last_block_len > 0 and chrom_name == ref["ref_name"]:
                logger.debug(f"TAIL {e} chrom_name {chrom_name} last_bin_pos {last_bin_pos:12,d} last_first_pos {last_first_pos:12,d} last_last_pos {last_last_pos:12,d} last_block_len {last_block_len:12,d}")
                extra_blocks.append({
                    "bin_n": -1,
                    "chunk_n": -1,
                    "real": lastReal,
                    "bytes": 0,
                    "block_len": last_block_len,
                    "bin_pos": last_bin_pos,
                    "first_pos": last_first_pos,
                    "last_pos": last_last_pos,
                })
            lastReal += last_block_len
            e += 1

        if len(extra_blocks) > 0:
            ref["last_block" ] = extra_blocks[-1]

        (n_intv,)  = get_values('<i')
        logger.debug(f"     n_intv              # 16kb intervals (for the linear index)         int32_t  {n_intv:15,d}")

        assert n_intv >     0, n_intv
        assert n_intv < 2**63, n_intv

        ref["n_intv"] = n_intv
        ioffs         = get_values('<' + ('Q'*n_intv))
        ioffs_be      = [None] * n_intv

        for intv_n in range(n_intv):
            # logger.debug(f"         intv_n {intv_n+1:15,d}/{n_intv:15,d}")
            ioff       = ioffs[intv_n]

            assert ioff >=     0
            assert ioff <  2**63

            ioff_real  = ioff >> 16
            ioff_bytes = ioff &  BLOCK_BYTES_MASK

            assert ioff_real  >=     0
            assert ioff_real  <  2**48

            assert ioff_bytes >=          0
            assert ioff_bytes <  BLOCK_SIZE
            
            ioff_bin_pos, ioff_first_pos, ioff_last_pos = -1, -1 ,-1
            if ioff in position_memoize:
                ioff_nfo = position_memoize[ioff]
                ioff_bin_pos, ioff_first_pos, ioff_last_pos, ioff_block_len = [ioff_nfo[k] for  k in ["bin_pos", "first_pos", "last_pos", "block_len"]]
            else:
                ioff_bin_pos, ioff_first_pos, ioff_last_pos, ioff_block_len, _, _, _, _ = getPos(inf, ioff_real, ioff_bytes, ref['ref_name'])
                ioff_nfo = {
                    "bin_n": -1,
                    "chunk_n": -1,
                    "real": ioff_real,
                    "bytes": ioff_bytes,
                    "block_len": ioff_block_len,
                    "bin_pos": ioff_bin_pos,
                    "first_pos": ioff_first_pos,
                    "last_pos": ioff_last_pos
                }
            
            ioffs_be[intv_n] = ioff_nfo

            if getLogLevel() == "DEBUG":
                logger.debug(f"======================== List of distinct intervals (n=n_intv [{n_intv:15,d}]) ================")
                logger.debug(f"         ioff            File offset of the first record in the interval uint64_t [{intv_n+1:15,d}/{n_intv:15,d}]")
                logger.debug(f"           real          File offset of the first record in the interval uint48_t {ioff_real:15,d}")
                logger.debug(f"           bytes         File offset of the first record in the interval uint16_t {ioff_bytes:15,d}")
                logger.debug(f"           block len                                                              {ioff_block_len:15,d}")
                logger.debug(f"           position 1st                                                           {ioff_first_pos:15,d}")
                logger.debug(f"           position bin                                                           {ioff_bin_pos:15,d}")
                logger.debug(f"           position lst                                                           {ioff_last_pos:15,d}")

        ref["intvs"] = ioffs_be

    n_no_coor = None
    try:
        (n_no_coor,)  = get_values('<Q')
    except Exception as e:
        # logger.debug(e)
        n_no_coor = None
        
    if n_no_coor is None:
        logger.debug(f" n_no_coor (optional)    # unmapped reads without coordinates set        uint64_t -")
    else:
        logger.debug(f" n_no_coor (optional)    # unmapped reads without coordinates set        uint64_t {n_no_coor:15,d}")

    data["n_no_coor"] = n_no_coor

    fhd.close()
    inf.close()

    logger.debug("finished reading")

    return data

if __name__ == "__main__":
    # logging.info("__main__")
    infile     = sys.argv[1]
    ingz, inid = get_filenames(infile)
    data       = read_tabix(ingz, compress=True, verbosity=0)

    setLogLevel(logger.debug)

    save(data, ingz)
