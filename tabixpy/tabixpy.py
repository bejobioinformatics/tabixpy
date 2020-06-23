import os
import sys
import gzip
import struct
import json
import logging

__format_ver__    = 4
__format_name__   = "TBJ"

COMPRESS          = True
BRICK_SIZE        = 216
BLOCK_SIZE        = 2**16
FILE_BYTES_MASK   = 0xFFFFFFFFFFFFF0000
BLOCK_BYTES_MASK  = 0x0000000000000FFFF


logging.basicConfig(
    level   = logging.DEBUG,
    format  = "%(asctime)8s - %(name)s - %(levelname)-6s - %(message)s",
    datefmt = "%H:%M:%S",
    # format="%(asctime)23s - %(name)s - %(levelname)6s - %(message)s",
    # datefmt='%Y-%m-%d %H:%M:%S',
)

logger        = logging.getLogger('tabixpy')
logger.setLevel(logging.INFO)

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
    else:
        ingz     = infile
        inid     = infile + ".tbi"

    return ingz, inid

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

MAX_BIN = (((1<<18)-1)//7)
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

def getPos(filehandle, real_pos, bytes_pos):
    # print(f"seek {real_pos} {bytes_pos}")
    filehandle.seek(real_pos, 0)

    bin_pos   = -1
    first_pos = -1
    last_pos  = -1
    with gzip.open(filehandle, 'rb') as g:
        block     = g.read(BLOCK_SIZE).decode()
        if len(block) > 0 and block[0] != "#":
            rows      = block.split("\n")
            first_row = None
            last_row  = None
            if len(rows) > 1:
                if len(rows[0]) == len(rows[1]):
                    first_row = rows[0]
                else:
                    first_row = rows[1]
                
                if len(rows[-2]) == len(rows[-1]):
                    last_row = rows[-1]
                else:
                    last_row = rows[-2]
            else:
                first_row = rows[0]
                last_row  = rows[-1]
            
            first_cols = first_row.split("\t")
            first_pos  = first_cols[1]
            first_pos  = int(first_pos)

            last_cols  = last_row.split("\t")
            last_pos   = last_cols[1]
            last_pos   = int(last_pos)

            bin_reg    = block[bytes_pos:bytes_pos+1024]
            bin_cols   = bin_reg.split("\t")
            bin_pos    = bin_cols[1]
            bin_pos    = int(bin_pos)
    return bin_pos, first_pos, last_pos


def read_tabix(infile):
    logger.info(f"reading {infile}")

    ingz, inid = get_filenames(infile)

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

        ref["n_bin"] = n_bin
        ref["bins" ] = [None] * n_bin

        logger.debug(f"======================== List of distinct bins (n=n_bin [{n_bin:15,d}]) ======================")
        last_chk_real_beg = None
        last_chk_real_end = None
        
        position_memoize = {}

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
                # "chunks" : [None] * n_chunk
            }
            # chunks = ref["bins"][bin_n]["chunks"]

            chunks = get_values('<' + ('Q'*(n_chunk*2)))

            assert all([i >     0 for i in chunks])
            assert all([i < 2**63 for i in chunks])

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

                assert chk_beg       < chk_end
                assert chk_real_beg  < chk_real_end
                # assert chk_bytes_beg < chk_bytes_end

                if last_chk_real_beg is not None:
                    assert last_chk_real_beg <  chk_real_beg
                    if last_chk_real_end is not None:
                        assert last_chk_real_end == chk_real_beg
                    
                if last_chk_real_end is not None:
                    assert last_chk_real_end <  chk_real_end
                    assert last_chk_real_end == chk_real_beg

                last_chk_real_beg = chk_real_beg
                last_chk_real_end = chk_real_end

                if chk_beg in position_memoize:
                    chk_bin_pos_beg, chk_first_pos_beg, chk_last_pos_beg = position_memoize[chk_beg]
                else:
                    chk_bin_pos_beg, chk_first_pos_beg, chk_last_pos_beg = getPos(inf, chk_real_beg, chk_bytes_beg)
                    position_memoize[chk_beg] = (chk_bin_pos_beg, chk_first_pos_beg, chk_last_pos_beg)

                if chk_end in position_memoize:
                    chk_bin_pos_end, chk_first_pos_end, chk_last_pos_end = position_memoize[chk_end]
                else:
                    chk_bin_pos_end, chk_first_pos_end, chk_last_pos_end = getPos(inf, chk_real_end, chk_bytes_end)
                    position_memoize[chk_end] = (chk_bin_pos_end, chk_first_pos_end, chk_last_pos_end)

                chunks_data["chunk_begin" ][chunk_n] = {
                    "real": chk_real_beg,
                    "bytes": chk_bytes_beg,
                    "bin_pos": chk_bin_pos_beg,
                    "first_pos": chk_first_pos_beg,
                    "last_pos": chk_last_pos_beg
                }
                
                chunks_data["chunk_end"   ][chunk_n] = {
                    "real": chk_real_end,
                    "bytes": chk_bytes_end,
                    "bin_pos": chk_bin_pos_end,
                    "first_pos": chk_first_pos_end,
                    "last_pos": chk_last_pos_end
                }

                if getLogLevel() == "DEBUG":
                    logger.debug(f"======================== List of chunks (n=n_chunk[{chunk_n:15,d}]) ============================")
                    logger.debug(f"             chunk_n {chunk_n+1:15,d}/{n_chunk:15,d}")
                    logger.debug(f"               cnk_beg   Virtual file offset of the start of the chunk   uint64_t {chk_real_beg:15,d} {chk_bytes_beg:15,d}")
                    logger.debug(f"               cnk_end   Virtual file offset of the end of the chunk     uint64_t {chk_real_end:15,d} {chk_bytes_end:15,d}")
                    logger.debug(f"               cnk_beg_1st_pos                                           uint64_t {chk_first_pos_beg:15,d}")
                    logger.debug(f"               cnk_end_1st_pos                                           uint64_t {chk_first_pos_end:15,d}")
                    logger.debug(f"               cnk_beg_bin_pos                                           uint64_t {chk_bin_pos_beg:15,d}")
                    logger.debug(f"               cnk_end_bin_pos                                           uint64_t {chk_bin_pos_end:15,d}")
                    logger.debug(f"               cnk_beg_lst_pos                                           uint64_t {chk_last_pos_beg:15,d}")
                    logger.debug(f"               cnk_end_lst_pos                                           uint64_t {chk_last_pos_end:15,d}")

            ref["bins"      ][bin_n]["chunks"] = chunks_data
            ref["bins_begin"]                  = chunks_data["chunk_begin"][ 0]
            ref["bins_end"  ]                  = chunks_data["chunk_begin"][-1]

        (n_intv,)  = get_values('<i')
        logger.debug(f"     n_intv              # 16kb intervals (for the linear index)         int32_t  {n_intv:15,d}")

        assert n_intv >     0, n_intv
        assert n_intv < 2**63, n_intv

        ref["n_intv"] = n_intv

        ioffs    = get_values('<' + ('Q'*n_intv))
        ioffs_be = [None] * n_intv

        for intv_n in range(n_intv):
            # logger.debug(f"         intv_n {intv_n+1:15,d}/{n_intv:15,d}")
            ioff       = ioffs[intv_n]

            assert ioff >     0
            assert ioff < 2**63

            ioff_real  = ioff >> 16
            ioff_bytes = ioff &  BLOCK_BYTES_MASK

            assert ioff_real  >=     0
            assert ioff_real  <  2**48

            assert ioff_bytes >=          0
            assert ioff_bytes <  BLOCK_SIZE
            
            ioff_bin_pos, ioff_first_pos, ioff_last_pos = -1, -1 ,-1
            if ioff in position_memoize:
                ioff_bin_pos, ioff_first_pos, ioff_last_pos = position_memoize[ioff]
            else:
                ioff_bin_pos, ioff_first_pos, ioff_last_pos = getPos(inf, ioff_real, ioff_bytes)
            
            ioffs_be[intv_n] = {
                "real": ioff_real,
                "bytes": ioff_bytes,
                "bin_pos": ioff_bin_pos,
                "first_pos": ioff_first_pos,
                "last_pos": ioff_last_pos
            }

            if getLogLevel() == "DEBUG":
                logger.debug(f"======================== List of distinct intervals (n=n_intv [{n_intv:15,d}]) ================")
                logger.debug(f"         ioff            File offset of the first record in the interval uint64_t [{intv_n+1:15,d}/{n_intv:15,d}]")
                logger.debug(f"           real          File offset of the first record in the interval uint48_t {ioff_real:15,d}")
                logger.debug(f"           bytes         File offset of the first record in the interval uint16_t {ioff_bytes:15,d}")
                logger.debug(f"           position 1st  File offset of the first record in the interval uint64_t {ioff_first_pos:15,d}")
                logger.debug(f"           position bin  File offset of the first record in the interval uint64_t {ioff_bin_pos:15,d}")
                logger.debug(f"           position lst  File offset of the first record in the interval uint64_t {ioff_last_pos:15,d}")

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
    infile     = sys.argv[1]
    ingz, inid = get_filenames(infile)
    data       = read_tabix(ingz, compress=True, verbosity=0)

    logger.setLevel(logging.DEBUG)

    save(data, ingz)
