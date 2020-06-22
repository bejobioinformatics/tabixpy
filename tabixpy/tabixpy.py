import os
import sys
import gzip
import struct
import json
import logging

__format_ver__    = 3
__format_name__   = "TBJ"

COMPRESS          = True

logging.basicConfig(
    level=logging.DEBUG,
    format="%(asctime)23s - %(name)s - %(levelname)6s - %(message)s"
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

def read_tabix(infile):
    logger.info(f"reading {infile}")

    ingz, inid = get_filenames(infile)

    assert os.path.exists(inid), inid

    fhd        = gzip.open(inid, "rb")
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

        (n_bin,)   = get_values('<i')
        
        logger.debug(f"     ref_n {ref_n+1:15,d}/{n_ref:15,d} ({names[ref_n]})")
        logger.debug(f"       n_bin             # distinct bins (for the binning index)         int32_t  {n_bin:15,d}")

        assert n_bin >     0, n_bin
        assert n_bin < 2**31, n_bin

        ref["n_bin"] = n_bin
        ref["bins" ] = [None] * n_bin

        logger.debug(f"======================== List of distinct bins (n=n_bin [{n_bin:15,d}]) ======================")
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

            block_bytes_mask    = 0x0000000000000FFFF
            real_file_offsets   = [c >> 16              for c in chunks]
            block_bytes_offsets = [c & block_bytes_mask for c in chunks]

            assert all([i >=     0 for i in real_file_offsets])
            assert all([i <  2**48 for i in real_file_offsets])

            assert all([i >=     0 for i in block_bytes_offsets])
            assert all([i <  2**16 for i in block_bytes_offsets])

            chunks_data = [None] * n_chunk

            for p, c in enumerate(range(0,len(chunks),2)):
                chk_beg       = chunks[c+0]
                chk_end       = chunks[c+1]

                chk_beg_real  = real_file_offsets[c+0]
                chk_end_real  = real_file_offsets[c+1]

                chk_beg_bytes = block_bytes_offsets[c+0]
                chk_end_bytes = block_bytes_offsets[c+1]

                # logger.debug(f"begin")
                # logger.debug(f"  block             {chk_beg:064b} {chk_beg:15,d}")
                # logger.debug(f"  mask              {block_bytes_mask:064b}")
                # logger.debug(f"  real offset       {chk_beg_real:064b} {chk_beg_real:15,d}")
                # logger.debug(f"  block byte offset {chk_beg_bytes:064b} {chk_beg_bytes:15,d}")
                # logger.debug(f"end")
                # logger.debug(f"  block             {chk_end:064b} {chk_end:15,d}")
                # logger.debug(f"  mask              {block_bytes_mask:064b}")
                # logger.debug(f"  real offset       {chk_end_real:064b} {chk_end_real:15,d}")
                # logger.debug(f"  block byte offset {chk_end_bytes:064b} {chk_end_bytes:15,d}")
                # logger.debug("")

                chunks_data[p] = [
                    [chk_beg_real, chk_beg_bytes],
                    [chk_end_real, chk_end_bytes]
                ]

            if getLogLevel() == "DEBUG":
                logger.debug(f"======================== List of chunks (n=n_chunk[{n_chunk:15,d}]) ============================")
                for chunk_n in range(n_chunk):
                    [
                        [chk_beg_real, chk_beg_bytes],
                        [chk_end_real, chk_end_bytes]
                    ] = chunks_data[chunk_n]
                    logger.debug(f"             chunk_n {chunk_n+1:15,d}/{n_chunk:15,d}")
                    logger.debug(f"               cnk_beg   Virtual file offset of the start of the chunk   uint64_t {chk_beg_real:15,d} {chk_beg_bytes:15,d}")
                    logger.debug(f"               cnk_end   Virtual file offset of the end of the chunk     uint64_t {chk_end_real:15,d} {chk_end_bytes:15,d}")
            ref["bins"][bin_n]["chunks"] = chunks_data

            # chunks = [(b,e) for b,e in zip(chunks[0::2], chunks[1::2])]
            # ref["bins"][bin_n]["chunks"] = chunks
            # if getLogLevel() == "DEBUG":
            #     logger.debug(f"======================== List of chunks (n=n_chunk[{n_chunk:15,d}]) ============================")
            #     for chunk_n in range(n_chunk):
            #         cnk_beg, cnk_end = chunks[chunk_n]
            #         logger.debug(f"             chunk_n {chunk_n+1:15,d}/{n_chunk:15,d}")
            #         logger.debug(f"               cnk_beg   Virtual file offset of the start of the chunk   uint64_t {cnk_beg:15,d}")
            #         logger.debug(f"               cnk_end   Virtual file offset of the end of the chunk     uint64_t {cnk_end:15,d}")

        (n_intv,)  = get_values('<i')
        logger.debug(f"     n_intv              # 16kb intervals (for the linear index)         int32_t  {n_intv:15,d}")

        assert n_intv >     0, n_intv
        assert n_intv < 2**63, n_intv

        ref["n_intv"] = n_intv

        ioffs = get_values('<' + ('Q'*n_intv))

        assert all([i >     0 for i in ioffs])
        assert all([i < 2**63 for i in ioffs])

        ref["intvs"] = ioffs

        if getLogLevel() == "DEBUG":
            logger.debug(f"======================== List of distinct intervals (n=n_intv [{n_intv:15,d}]) ================")
            for intv_n in range(n_intv):
                # logger.debug(f"         intv_n {intv_n+1:15,d}/{n_intv:15,d}")
                ioff = ioffs[intv_n]
                logger.debug(f"         ioff            File offset of the first record in the interval uint64_t {ioff:15,d} [{intv_n+1:15,d}/{n_intv:15,d}]")

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

    logger.debug("finished reading")

    return data

if __name__ == "__main__":
    infile     = sys.argv[1]
    ingz, inid = get_filenames(infile)
    data       = read_tabix(ingz, compress=True, verbosity=0)

    logger.setLevel(logging.DEBUG)

    save(data, ingz)
