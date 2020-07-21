import os
import struct
import json

from ._gzip       import gzip, getBlock, EOF
from ._io         import getFilenames, genStructValueGetter
from ._logger     import logger, getLogLevel
from ._consts     import (
    TABIX_EXTENSION,
    TABIX_MAGIC,
    TABIX_FILE_BYTES_MASK,
    TABIX_BLOCK_BYTES_MASK,
    TABIX_MAX_BIN
)

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
    #define TABIX_MAX_BIN (((1<<18)-1)/7)
    
    res = [None] * TABIX_MAX_BIN
    
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

def parseBlock(block, bytes_pos, chrom):
    bin_pos        = -1
    first_pos      = -1
    last_pos       = -1
    first_chrom    = None
    len_first_cols = -1
    num_rows       = -1

    assert len(block) > 0, f"empty block '{block}'"
    assert bytes_pos  < len(block)

    first_row = None
    last_row  = None
    rows      = [row.split("\t") for row in block.split("\n") if len(row) > 0 and row[0] != "#"]
    num_rows  = len(rows)

    if num_rows == 0:
        if getLogLevel() == "DEBUG":
            logger.debug(f"num_rows {num_rows} == 0")
        return -1, -1, -1, None, -1, -1

    if num_rows < 3:
        if getLogLevel() == "DEBUG":
            logger.debug(f"num_rows {num_rows} < 3")
        return -1, -1, -1, None, -1, -1

    rows       = [row for row in rows if len(row) == len(rows[1])]

    if num_rows == 0:
        if getLogLevel() == "DEBUG":
            logger.debug(f"num_rows {num_rows} == 0")
        return -1, -1, -1, None, -1, -1

    if num_rows < 3:
        if getLogLevel() == "DEBUG":
            logger.debug(f"num_rows {num_rows} < 3")
        return -1, -1, -1, None, -1, -1
    
    # rows_group = list(set([row[0] for row in rows]))
    # rows_data  = []
    # logger.info(rows_group)

    if chrom is not None:
        rows      = [row for row in rows if row[0] == chrom]

    num_rows  = len(rows)

    if num_rows == 0:
        if getLogLevel() == "DEBUG":
            logger.debug(f"num_rows {num_rows} == 0")
        return -1, -1, -1, None, -1, -1

    if num_rows < 3:
        if getLogLevel() == "DEBUG":
            logger.debug(f"num_rows {num_rows} < 3")
        return -1, -1, -1, None, -1, -1

    elif num_rows == 1:
        first_row = rows[0]
        last_row  = rows[0]
        
        if getLogLevel() == "DEBUG":
            logger.debug(f"num_rows {num_rows} == 1 :: first_row {first_row[:2]} last_row {last_row[:2]}")
        
        first_cols = first_row
        last_cols  = last_row

    else:
        first_cols  = rows[ 0]
        second_cols = rows[ 1]

        if len(first_cols) != len(second_cols) or first_cols[0] != second_cols[0]:
            first_cols  = second_cols
            num_rows   -= 1
        
        assert chrom is None or first_cols[0] == chrom
        
        for i in range(len(rows)):
            if getLogLevel() == "DEBUG":
                logger.debug(f"i {i} (i*-1)-1 {(i*-1)-1:3d} len(rows) {len(rows):3d}")
            
            last_cols    = rows[(i*-1)-1]

            if len(last_cols) >= 2 and last_cols[0] == first_cols[0] and len(last_cols[1]) > 0:
                try:
                    last_cols_v  = int(last_cols[1])
                except ValueError as e:
                    raise ValueError(f"invalid last_col position {last_cols} :: {e}")
                
                try:
                    first_cols_v = int(first_cols[1])
                except ValueError as e:
                    raise ValueError(f"invalid fist_col position {first_cols} :: {e}")
            
                if last_cols_v > first_cols_v:
                    if chrom is None or last_cols[0] == chrom: 
                        # if getLogLevel() == "DEBUG":
                        #     logger.debug(f"i {i} (i*-1)-1 {(i*-1)-1:3d} len(rows) {len(rows):3d} num_rows {num_rows:3d}\n\tlen(last_cols) {len(last_cols):3d} == len(first_cols) {len(first_cols):3d} or last_cols[0] {last_cols[0] or None} == first_cols[0] {first_cols[0] or None} - last_cols[1] {last_cols[1] or None} == first_cols[1] {first_cols[1] or None}\n")
                        break
            else:
                # if getLogLevel() == "DEBUG":
                #     logger.debug(f"i {i} (i*-1)-1 {(i*-1)-1:3d} len(rows) {len(rows):3d} num_rows {num_rows:3d}\n\tlen(last_cols) {len(last_cols):3d} != len(first_cols) {len(first_cols):3d} or last_cols {last_cols} first_cols {first_cols}")
                num_rows   -= 1

    len_first_cols = len(first_cols)
    len_last_cols  = len(last_cols)

    if len_first_cols < 2:
        return -1, -1, -1, None, -1, -1

    if len_first_cols < 2:
        return -1, -1, -1, None, -1, -1

    if len_first_cols < len_last_cols:
        return -1, -1, -1, None, -1, -1

    assert len_first_cols >= len_last_cols, f"{len_first_cols} >= {len_last_cols}\n{len_first_cols} {first_cols}\n{len_last_cols} {last_cols}"
    assert len_first_cols >= 2            , f"{len_first_cols} >= 2\n{len_first_cols} {first_cols}"
    assert len_last_cols  >= 2            , f"{len_last_cols } >= 2\n{len_last_cols}  {last_cols }"

    first_chrom = first_cols[0]
    last_chrom  = last_cols[0]

    assert chrom is None or first_chrom == last_chrom, f"first_chrom {first_chrom} == last_chrom {last_chrom}"

    first_pos   = first_cols[1]
    try:
        first_pos   = int(first_pos)
    except:
        raise ValueError(f"invalid positions {first_pos} :: {first_row}")

    last_pos    = last_cols[1]
    try:
        last_pos    = int(last_pos)
    except:
        raise ValueError(f"invalid positions {last_pos} :: {last_row}")

    bin_reg      = block[bytes_pos:].split("\n")[0]
    bin_cols     = bin_reg.split("\t")
    len_bin_cols = len(bin_cols)
    if len_bin_cols == len_first_cols:
        bin_chrom  = bin_cols[0]
        bin_pos    = bin_cols[1]
        
        assert len_bin_cols == len_first_cols, f"len(bin_cols) {len_bin_cols} == len(first_cols) {len_first_cols} :: bytes_pos {bytes_pos} chrom {chrom}\nblock[         :100] {block[:100]}\nblock[bytes_pos:   ] {block[bytes_pos:bytes_pos+100]}"

        try:
            bin_pos    = int(bin_pos)
        except ValueError as e:
            logger.error(e)
            logger.error(bin_reg)
            raise

        assert len_bin_cols >= len_last_cols, f"{len_bin_cols} == {len_last_cols}\n{len_bin_cols} {bin_cols}\n{len_last_cols} {last_cols}"

    return bin_pos, first_pos, last_pos, first_chrom, len_first_cols, num_rows

def getPos(filehandle, real_pos, bytes_pos, chrom):
    # logger.debug(f"getPos :: real_pos {real_pos} bytes_pos {bytes_pos} chrom {chrom}")

    if filehandle is None:
        raise ValueError("Filehandle is empty")

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

def getPosData(bin_n, chunk_n, inf, real_pos, bytes_pos, chrom):
    if inf is not None:
        bin_pos, first_pos, last_pos, block_len, _, _, _, _ = getPos(inf, real_pos, bytes_pos, chrom)

        chunk_nfo = {
            "bin_n": bin_n,
            "chunk_n": chunk_n,
            "real": real_pos,
            "bytes": bytes_pos,
            "block_len": block_len,
            "bin_pos": bin_pos,
            "first_pos": first_pos,
            "last_pos": last_pos
        }
    else:
        chunk_nfo = {
            "bin_n": bin_n,
            "chunk_n": chunk_n,
            "real": real_pos,
            "bytes": bytes_pos
        }


    return chunk_nfo

def readTabix(infile):
    logger.info(f"reading {infile}{TABIX_EXTENSION}")

    (ingz, inid, inbj, inbk) = getFilenames(infile)

    assert os.path.exists(inid), inid

    fhd        = gzip.open(inid, "rb")

    inf = None
    if os.path.exists(ingz):
        inf        = open(ingz, "rb")
    
    get_values = genStructValueGetter(fhd)
    data       = {}

    (magic,)   = get_values('<4s')

    assert magic == TABIX_MAGIC, f'invalid magic: {magic}'

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

            # real_file_offsets   = [c & TABIX_FILE_BYTES_MASK   for c in chunks]
            real_file_offsets   = [c >> 16              for c in chunks]
            block_bytes_offsets = [c & TABIX_BLOCK_BYTES_MASK for c in chunks]

            assert all([i >=     0 for i in real_file_offsets])
            assert all([i <  2**48 for i in real_file_offsets])

            assert all([i >=     0 for i in block_bytes_offsets])
            assert all([i <  2**16 for i in block_bytes_offsets])

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
                # logger.debug(f"  mask              {TABIX_BLOCK_BYTES_MASK:064b}")
                # logger.debug(f"  real offset       {chk_real_beg:064b} {chk_real_beg:15,d}")
                # logger.debug(f"  block byte offset {chk_bytes_beg:064b} {chk_bytes_beg:15,d}")
                # logger.debug(f"end")
                # logger.debug(f"  block             {chk_end:064b} {chk_end:15,d}")
                # logger.debug(f"  mask              {TABIX_BLOCK_BYTES_MASK:064b}")
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

                chunk_nfo_beg = position_memoize.get(chk_beg, None)
                if chunk_nfo_beg is None:
                    chunk_nfo_beg = getPosData(bin_n, chunk_n, inf, chk_real_beg, chk_bytes_beg, ref['ref_name'])
                    position_memoize[chk_beg] = chunk_nfo_beg

                chk_bin_pos_beg, chk_first_pos_beg, chk_last_pos_beg, chk_block_len_beg = [chunk_nfo_beg.get(k, None) for  k in ["bin_pos", "first_pos", "last_pos", "block_len"]]


                chunk_nfo_end = position_memoize.get(chk_end, None)
                if chunk_nfo_end is None:
                    chunk_nfo_end = getPosData(bin_n, chunk_n, inf, chk_real_end, chk_bytes_end, ref['ref_name'])
                    position_memoize[chk_end] = chunk_nfo_end

                chk_bin_pos_end, chk_first_pos_end, chk_last_pos_end, chk_block_len_end = [chunk_nfo_end.get(k, None) for  k in ["bin_pos", "first_pos", "last_pos", "block_len"]]


                chunks_data["chunk_begin" ][chunk_n] = chunk_nfo_beg
                chunks_data["chunk_end"   ][chunk_n] = chunk_nfo_end

                if inf is None:
                    if chunk_nfo_beg["real"] < bin_min_pos:
                        ref["bins_begin"] = chunks_data["chunk_begin" ][chunk_n]
                        bin_min_pos = chunk_nfo_beg["real"]

                    if chunk_nfo_end["real"] > bin_max_pos:
                        ref["bins_end"  ] = chunks_data["chunk_end"   ][chunk_n]
                        bin_max_pos = chunk_nfo_end["real"]

                else:
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
        last_bin_pos, last_first_pos, last_last_pos, last_block_len = [last_bin.get(k, None) for k in ['bin_pos', 'first_pos', 'last_pos', 'block_len']]
        logger.debug(f"LAST 0 last_bin_pos {last_bin_pos or -1:12,d} last_first_pos {last_first_pos or -1:12,d} last_last_pos {last_last_pos or -1:12,d} last_block_len {last_block_len or -1:12,d}")

        ref["first_block"] = ref["bins_begin"]
        ref["last_block" ] = ref["bins_end"  ]

        lastReal = None
        if inf is not None:
            lastReal     = last_bin["real"] + last_block_len
        
            chrom_name   = ref["ref_name"]
            extra_blocks = []
            e            = 0

            while last_block_len > 0 and last_last_pos > 0 and chrom_name == ref["ref_name"]:
                chrom_name = ref['ref_name']
                # last_bin_pos, last_first_pos, last_last_pos, last_block_len, _, _, _, _ = getPos(inf, lastReal, 0, ref['ref_name'])
                last_nfo       = getPosData(-1, -1, inf, lastReal, 0, ref['ref_name'])

                last_bin_pos, last_first_pos, last_last_pos, last_block_len = [last_nfo[k] for  k in ["bin_pos", "first_pos", "last_pos", "block_len"]]

                if last_block_len > 0 and last_last_pos > 0:# and chrom_name == ref["ref_name"]:
                    logger.debug(f"TAIL {e} chrom_name {chrom_name} last_bin_pos {last_bin_pos:12,d} last_first_pos {last_first_pos:12,d} last_last_pos {last_last_pos:12,d} last_block_len {last_block_len:12,d}")
                    extra_blocks.append(last_nfo)
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
            ioff_bytes = ioff &  TABIX_BLOCK_BYTES_MASK

            assert ioff_real  >=     0
            assert ioff_real  <  2**48

            assert ioff_bytes >=     0
            assert ioff_bytes <  2**16
            
            ioff_nfo = position_memoize.get(ioff, None)
            if ioff_nfo is None:
                ioff_nfo = getPosData(-1 ,-1, inf, ioff_real, ioff_bytes, ref['ref_name'])

            ioff_bin_pos, ioff_first_pos, ioff_last_pos, ioff_block_len = [ioff_nfo.get(k, None) for  k in ["bin_pos", "first_pos", "last_pos", "block_len"]]
            
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
    if inf:
        inf.close()

    logger.debug("finished reading")

    return data
