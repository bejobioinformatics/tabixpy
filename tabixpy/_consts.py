BLOCK_SIZE          = 2**16
GZIP_MAGIC          = b"\x1F\x8B"
COMPRESS            = True

TABIX_FORMAT_NAME   = "TBI"
TABIX_EXTENSION     = '.tbi'
TABIX_MAGIC         = b'TBI\x01'

TABIXPY_FORMAT_VER  = 5
TABIXPY_FORMAT_NAME = "TBJ"
TABIXPY_EXTENSION   = '.tbj'

VCFBGZ_FORMAT_VER   = 1
VCFBGZ_FORMAT_NAME  = "TBK"
VCFBGZ_EXTENSION    = ".tbk"
VCFBGZ_EOF          = bytes.fromhex('000102030405060708090A0B0C0D0E0F')

# samtools format specs:
# https://samtools.github.io/hts-specs/SAMv1.pdf
# https://github.com/xbrianh/bgzip/blob/master/bgzip/__init__.py
BGZIP_EOF = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")
