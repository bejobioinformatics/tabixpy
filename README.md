tabixpy
=======

Tabix parser writtern in Python3.

CI
--
![Upload Python Package](https://github.com/bejobioinformatics/tabixpy/workflows/Upload%20Python%20Package/badge.svg)

Install
-------

pip install tabixpy

Tabix
-----

https://samtools.github.io/hts-specs/tabix.pdf


```
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
======================= List of indices (n=n_ref )            =======================
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
```

Offset calculation

```python
block_bytes_mask    = 0x0000000000000FFFF
real_file_offset    = chunk >> 16
block_bytes_offset  = chunk & block_bytes_mask
```

```
begin
    block             0000000000000000000000000110100011010101100111111001011010001111 450,260,604,559
    mask              0000000000000000000000000000000000000000000000001111111111111111
    real offset       0000000000000000000000000000000000000000011010001101010110011111       6,870,431
    block byte offset 0000000000000000000000000000000000000000000000001001011010001111          38,543
end
    block             0000000000000000000000000110101001010101110101110011001010010100 456,706,699,924
    mask              0000000000000000000000000000000000000000000000001111111111111111
    real offset       0000000000000000000000000000000000000000011010100101010111010111       6,968,791
    block byte offset 0000000000000000000000000000000000000000000000000011001010010100          12,948
```

BGZIP
-----

http://samtools.github.io/hts-specs/SAMv1.pdf


```
The random access method to be described next limits the uncompressed contents
of each BGZF block to a maximum of 216 bytes of data. Thus while ISIZE is
stored as a uint32 t as per the gzip format, in BGZF it is limited to the range
[0, 65536]. BSIZE can represent BGZF block sizes in the range [1, 65536],
though typically BSIZE will be rather less than ISIZE due to compression.

4.1.1 Random access
BGZF files support random access through the BAM file index. To achieve this,
the BAM file index uses virtual file offsets into the BGZF file. Each virtual
file offset is an unsigned 64-bit integer, defined as: coffset<<16|uoffset,
where coffset is an unsigned byte offset into the BGZF file to the beginning of
a BGZF block, and uoffset is an unsigned byte offset into the uncompressed data
stream represented by that BGZF block. Virtual file offsets can be compared,
but subtraction between virtual file offsets and addition between a virtual
offset and an integer are both disallowed.
```

TABIX
-----

https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3042176/
https://samtools.github.io/hts-specs/tabix.pdf

```
2.1 Sorting and BGZF compression
Bgzip compresses the data file in the BGZF format, which is the concatenation
of a series of gzip blocks with each block holding at most 216 bytes of
uncompressed data.

In the compressed file, each uncompressed byte in the text data file is
assigned a unique 64-bit virtual file offset where the higher 48 bits keep the
real file offset of the start of the gzip block the byte falls in, and the
lower 16 bits store the offset of the byte inside the gzip block. Given a
virtual file offset, one can directly seek to the start of the gzip block using
the higher 48 bits, decompress the block and retrieve the byte pointed by the
lower 16 bits of the virtual offset. Random access can thus be achieved without
the help of additional index structures. As gzip works with concatenated gzip
files, it can also seamlessly decompress a BGZF file. The detailed description
of the BGZF format is described in the SAM specification.

2.2 Coupled binning and linear indices
Tabix builds two types of indices for a data file: a binning index and a linear
index. We can actually achieve fast retrieval with only one of them. However,
using the binning index alone may incur many unnecessary seek calls, while
using the linear index alone has bad worst-case performance (when some records
span very long distances). Using them together avoids their weakness.

2.2.1 The binning index
The basic idea of binning is to cluster records into large intervals, called
bins. A record is assigned to bin k if k is the bin of the smallest size that
fully contains the record. For each bin, we keep in the index file the virtual
file offsets of all records assigned to the bin. When we search for records
overlapping a query interval, we first collect bins overlapping the interval
and then test each record in the collected bins for overlaps.

In principle, bins can be selected freely as long as each record can be
assigned to a bin. In the Tabix binning index, we adopt a multilevel binning
scheme where bins at same level are non-overlapping and of the same size. In
Tabix, each bin k, 0<=k<=37,449, represents a half-close-half-open interval
[(k-ol)sl, (k-ol+1)sl), where l = [log2(7k+1)/3] is the level of the bin,
sl = 2^(29-3l) is the size of the bin at level l and ol = (2^3l - 1)/7 is the
offset at l. In this scheme, bin 0 spans 512 Mb, 1-8 span 64 Mb, 9-72 8 Mb,
73-584 1 Mb, 585-4680 128 kb and 4681-37449 span 16 kb intervals. The scheme is
very similar to the UCSC binning (Kent et al., 2002) except that in UCSC,
0<=k<=4681 and therefore the smallest bin size is 128 kb.

2.2.2 The linear index
With the binning index alone, we frequently need to seek to records assigned to
bins at lower levels, in particular bin 0, the bin spanning the entire
sequence. However, frequently records at lower levels do not overlap the query
interval, which leads to unsuccessful seeks and hurts performance especially
when data are transferred over network. The linear index is to reduce
unnecessary seek calls in this case (Table 1).

In the linear index, we keep for each tiling 16 kb window the virtual file
offset of the leftmost record (i.e. having the smallest start coordinate) that
overlaps the window. When we search for records overlapping a query interval,
we will know from the index the leftmost record that possibly overlaps the
query interval. Records having smaller coordinates than this leftmost record
can be skipped and unsuccessful seek calls can be saved.

When the data files are sorted by coordinates, records assigned to the same bin
tend to be adjacent. Thus in the index file, we do not need to keep the virtual
file offset of each record, but only to keep the start offset of a chunk of
records assigned to the same bin.
```

Schema
------

https://jsonschema.net/home
https://www.jsonschemavalidator.net/


Example output
--------------

JSON

```JSON
{
 "n_ref": 1,
 "format": 2,
 "col_seq": 1,
 "col_beg": 2,
 "col_end": 0,
 "meta": "#",
 "skip": 0,
 "l_nm": 11,
 "names": [
  "SL2.50ch00"
 ],
 "refs": [
  {
   "ref_n": 0,
   "ref_name": "SL2.50ch00",
   "n_bin": 86,
   "bins": [
    {
     "bin_n": 0,
     "bin": 4681,
     "n_chunk": 1,
     "chunks": {
      "chunk_begin": [
       {
        "real": 0,
        "bytes": 29542,
        "bin_pos": -1,
        "first_pos": -1,
        "last_pos": -1
       }
      ],
      "chunk_end": [
       {
        "real": 124525,
        "bytes": 19630,
        "bin_pos": 16388,
        "first_pos": 16141,
        "last_pos": 17808
       }
      ]
     }
    }
    ],
   "bins_begin": {
    "real": 7021611,
    "bytes": 4631,
    "bin_pos": 1392700,
    "first_pos": 1392519,
    "last_pos": 1393971
   },
   "bins_end": {
    "real": 7021611,
    "bytes": 4631,
    "bin_pos": 1392700,
    "first_pos": 1392519,
    "last_pos": 1393971
   },
   "n_intv": 86,
   "intvs": [
    {
     "real": 7021611,
     "bytes": 4631,
     "bin_pos": 1392700,
     "first_pos": 1392519,
     "last_pos": 1393971
    }
   ]
  }
 ],
 "n_no_coor": null,
 "__format_name__": "TBJ",
 "__format_ver__": 4
}
```

Timming
-------

```
2020-06-12 11:30:34,716 - tabixpy -   INFO - reading annotated_tomato_150.100000.vcf.gz
2020-06-12 11:30:34,738 - tabixpy -   INFO - saving  annotated_tomato_150.100000.vcf.gz.tbj
                   ,024
2020-06-12 11:31:16,506 - tabixpy -   INFO - reading annotated_tomato_150.vcf.bgz
2020-06-12 11:31:24,152 - tabixpy -   INFO - saving  annotated_tomato_150.vcf.bgz.tbj
                  8,646
```

File Sizes
----------

Compressed

```
1.1K annotated_tomato_150.100000.vcf.gz.tbi
2.0K annotated_tomato_150.100000.vcf.gz.tbj
727K annotated_tomato_150.vcf.bgz.tbi
1.2M annotated_tomato_150.vcf.bgz.tbj
```

Uncompressed

```
1.1K annotated_tomato_150.100000.vcf.gz.tbi
 15K annotated_tomato_150.100000.vcf.gz.tbj
727K annotated_tomato_150.vcf.bgz.tbi
8.4M annotated_tomato_150.vcf.bgz.tbj
```
