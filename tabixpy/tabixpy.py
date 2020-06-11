import os
import sys
import gzip
import struct
import json
import logging

__format_ver__    = 2
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

            chunks = [(b,e) for b,e in zip(chunks[0::2], chunks[1::2])]

            ref["bins"][bin_n]["chunks"] = chunks
            
            if getLogLevel() == "DEBUG":
                logger.debug(f"======================== List of chunks (n=n_chunk[{n_chunk:15,d}]) ============================")
                for chunk_n in range(n_chunk):
                    cnk_beg, cnk_end = chunks[chunk_n]
                    logger.debug(f"             chunk_n {chunk_n+1:15,d}/{n_chunk:15,d}")
                    logger.debug(f"               cnk_beg   Virtual file offset of the start of the chunk   uint64_t {cnk_beg:15,d}")
                    logger.debug(f"               cnk_end   Virtual file offset of the end of the chunk     uint64_t {cnk_end:15,d}")

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


"""
https://samtools.github.io/hts-specs/tabix.pdf

Field                   Description                                     Type     Value
---------------------------------------------------------------------------------------
magic                   Magic string                                    char[4]  TBI\1
n_ref                   # sequences                                     int32_t
format                  Format (0: generic; 1: SAM; 2: VCF)             int32_t
col_seq                 Column for the sequence name                    int32_t
col_beg                 Column for the start of a region                int32_t
col_end                 Column for the end of a region                  int32_t
meta                    Leading character for comment lines             int32_t
skip                    # lines to skip at the beginning                int32_t
l_nm                    Length of concatenated sequence names           int32_t
names                   Concatenated names, each zero terminated        char[l_nm]
======================= List of indices (n=n_ref)             =======================
    n_bin               # distinct bins (for the binning index)         int32_t
======================= List of distinct bins (n=n_bin)       =======================
        bin             Distinct bin number                             uint32_t
        n_chunk         # chunks                                        int32_t
======================= List of chunks (n=n_chunk)            =======================
            cnk_beg     Virtual file offset of the start of the chunk   uint64_t
            cnk_end     Virtual file offset of the end of the chunk     uint64_t
    n_intv              # 16kb intervals (for the linear index)         int32_t
======================= List of distinct intervals (n=n_intv) =======================
        ioff            File offset of the first record in the interval uint64_t
n_no_coor (optional)    # unmapped reads without coordinates set        uint64_t

Notes:
- The index file is BGZF compressed.

- All integers are little-endian.

- When (format&0x10000) is true, the coordinate follows the BED rule (i.e. half-closed-half-open and
zero based); otherwise, the coordinate follows the GFF rule (closed and one based).

- For the SAM format, the end of a region equals POS plus the reference length in the alignment, inferred
from CIGAR. For the VCF format, the end of a region equals POS plus the size of the deletion.

- Field col beg may equal col end, and in this case, the end of a region is end=beg+1.

- Example:
  For GFF, format=0      , col seq=1, col beg=4, col end=5, meta=‘#’ and skip=0.
  For BED, format=0x10000, col seq=1, col beg=2, col end=3, meta=‘#’ and skip=0.

- Given a zero-based, half-closed and half-open region [beg, end), the bin number is calculated with
the following C function:
    int reg2bin(int beg, int end) {
        --end;
        if (beg>>14 == end>>14) return ((1<<15)-1)/7 + (beg>>14);
        if (beg>>17 == end>>17) return ((1<<12)-1)/7 + (beg>>17);
        if (beg>>20 == end>>20) return ((1<< 9)-1)/7 + (beg>>20);
        if (beg>>23 == end>>23) return ((1<< 6)-1)/7 + (beg>>23);
        if (beg>>26 == end>>26) return ((1<< 3)-1)/7 + (beg>>26);
        return 0;
    }

- The list of bins that may overlap a region [beg, end) can be obtained with the following C function:
    #define MAX_BIN (((1<<18)-1)/7)
    int reg2bins(int rbeg, int rend, uint16_t list[MAX_BIN]) {
        int i = 0, k;
        --rend;
        list[i++] = 0;
        for (k =    1 + (rbeg>>26); k <=    1 + (rend>>26); ++k) list[i++] = k;
        for (k =    9 + (rbeg>>23); k <=    9 + (rend>>23); ++k) list[i++] = k;
        for (k =   73 + (rbeg>>20); k <=   73 + (rend>>20); ++k) list[i++] = k;
        for (k =  585 + (rbeg>>17); k <=  585 + (rend>>17); ++k) list[i++] = k;
        for (k = 4681 + (rbeg>>14); k <= 4681 + (rend>>14); ++k) list[i++] = k;
        return i; // #elements in list[]
    }
"""
