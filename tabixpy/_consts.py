BLOCK_SIZE        = 2**16
GZIP_MAGIC        = b"\x1F\x8B"
COMPRESS          = True

TABIXPY_FORMAT_VER  = 5
TABIXPY_FORMAT_NAME = "TBJ"


# samtools format specs:
# https://samtools.github.io/hts-specs/SAMv1.pdf
# https://github.com/xbrianh/bgzip/blob/master/bgzip/__init__.py
BGZIP_EOF = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")
