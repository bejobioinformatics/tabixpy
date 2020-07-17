import gzip
import struct

from ._logger     import logger, getLogLevel
from ._consts     import BLOCK_SIZE

GZIP_MAGIC          = b"\x1F\x8B"

# samtools format specs:
# https://samtools.github.io/hts-specs/SAMv1.pdf
# https://github.com/xbrianh/bgzip/blob/master/bgzip/__init__.py
BGZIP_EOF = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")


class EOF:
    __slot__ = []
    pass


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

            if getLogLevel() == "DEBUG":
                logger.debug(f"bn {bn} len(lines) {len(lines)} :: filter")
            lines   = [line for line in lines if len(line) > 0 and line[0] != "#"] # filter out empty and comment lines 
            columns = [line.split("\t") for line in lines]

            if len(columns[-1]) != len(columns[1]): #incomplete last line
                if getLogLevel() == "DEBUG":
                    logger.debug(f"bn {bn} len(lines) {len(lines)} :: pop")
                lastLine = lines.pop()
                columns.pop()

            if len(columns[0]) != len(columns[1]): #incomplete first line
                if getLogLevel() == "DEBUG":
                    logger.debug(f"bn {bn} len(lines) {len(lines)} :: shift")
                lines   = lines[1:]
                columns = columns[1:]

            if numCols is None:
                numCols = len(columns[0])

            if self._chrom is not None:
                if getLogLevel() == "DEBUG":
                    logger.debug(f"bn {bn} len(lines) {len(lines)} :: chrom {self._chrom}")
                lines   = [line for ln, line in enumerate(lines) if columns[ln][0] == self._chrom]
                columns = [col  for col      in columns          if col[0]         == self._chrom]

            if len(lines) == 0:
                if getLogLevel() == "DEBUG":
                    logger.debug(f"bn {bn} len(lines) == 0 :: return")
                return

            for ln, line in enumerate(lines):
                cols = columns[ln]
                if getLogLevel() == "DEBUG":
                    logger.debug(f"bn {bn} len(lines) == {len(lines)} :: chrom {self._chrom} begin {self._begin} end {self._end} :: ln {ln}")

                if self._begin is not None or self._end is not None:
                    pos = int(cols[1])
                    
                    if self._begin is not None and pos < self._begin:
                        if getLogLevel() == "DEBUG":
                            logger.debug(f" pos {pos} < self._begin {self._begin}")
                        continue

                    if self._end is not None and pos >= self._end:
                        if getLogLevel() == "DEBUG":
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



def readGzipHeader(fp):
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

def getBlock(filehandle, real_pos, block_len=None):
    filehandle.seek(real_pos, 0)

    if block_len is None:
        block_len = readGzipHeader(filehandle)

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
    assert len(block)     == len(blocktext)

    # logging.error(f"len(block) {len(blockdata)} == header_len {header_len} + len(bgzip_eof) {len(bgzip_eof)}")
    if len(blockdata) == len(BGZIP_EOF):
        if blockdata == BGZIP_EOF:
            return EOF(), -1

    # logger.debug(f" header_len  {header_len:12,d}")
    # logger.debug(f" blockdata   {len(blockdata):12,d}")
    # logger.debug(f" block       {len(block):12,d}")
    # logger.debug(f" blocktext   {len(blocktext):12,d}")

    assert len(blockdata) == block_len, f"block data has wrong size :: real_pos {real_pos} block_len {block_len}"
    assert len(blockdata) > 0         , f"block data has size zero :: real_pos {real_pos} block_len {block_len}"
    assert len(block)     > 0         , f"block has size zero :: real_pos {real_pos} block_len {block_len}"
    assert len(blocktext) > 0         , f"block text has size zero :: real_pos {real_pos} block_len {block_len}"

    # rows = block.decode().split("\n")

    # logging.info(len(rows))
    # logging.info(rows[0])
    # logging.info(rows[-2])
    # logging.info(rows[-1])

    return blocktext, block_len
