import os
import sys

DEBUG     = False
COMPRESS  = True
OVERWRITE = True

sys.path.insert(0, '../..')

import tabixpy

def runTest(testName, infile, expects, indexType):
    tb         = tabixpy.Tabix(infile, indexType=indexType)

    gzfile     = tb.bgz
    indexFile  = tb.indexFile
    sourceFile = tb.sourceFile

    if DEBUG:
        tabixpy.setLogLevel(tabixpy.logging.INFO)

    if os.path.exists(indexFile):
        tabixpy.logger.info(f"index file {indexFile} already exists")
        
    else:
        tabixpy.logger.info(f"index file {indexFile} does not exists. creating")

        tabix1 = tabixpy.Tabix(gzfile, indexType=indexType)
        tabix1.load()
        # tabix1.save(compress=COMPRESS, overwrite=OVERWRITE)

        tabixpy.logger.info(f"index file {indexFile} does not exists. reloading")

        tabix2 = tabixpy.Tabix(gzfile, indexType=indexType)
        tabix2.load()

        assert tabix1.data == tabix2.data

        del tabix1
        del tabix2

        tabixpy.logger.info(f"index file {indexFile} does not exists. removing")

    if DEBUG:
        tabixpy.setLogLevel(tabixpy.logging.DEBUG)

    if not os.path.exists(gzfile):
        #TODO: implement
        tabixpy.logger.warning(f"source file {gzfile} does not exists")
        return

    tabixpy.logger.info(f"loading index file {indexFile}")
    tb.load()

    tabixpy.logger.debug(tb.chromosomes)

    for test_num, (chrom_idx, begin, end, minPos, maxPos, count) in enumerate(expects):
        chrom = tb.chromosomes[chrom_idx]

        tabixpy.logger.info(f"testName {testName} test_num {test_num} chrom_idx {chrom_idx} chrom {chrom}, begin {begin}, end {end}, minPos {minPos}, maxPos {maxPos}")
        
        vals = list(tb.getChromosomeIter(chrom, begin=begin, end=end))

        if len(vals) > 0:
            assert int(vals[ 0][1]) == minPos, f"{vals[ 0][1]} == {minPos}"
            assert int(vals[-1][1]) == maxPos, f"{vals[-1][1]} == {maxPos}"
            if len(vals) > 10:
                for row_num, row in enumerate(vals[:5]):
                    tabixpy.logger.info(f"{row_num+1:6d} {row[:2]}")
                for row_num, row in enumerate(vals[-5:]):
                    tabixpy.logger.info(f"{row_num+1+len(vals)-5:6d} {row[:2]}")
            else:
                for row_num, row in enumerate(vals):
                    tabixpy.logger.info(f"{row_num+1:6d} {row[:2]}")

        assert len(vals) == count, f"{len(vals)} == {count}"

def runTests(tests, indexTypes):
    for indexType in indexTypes:
        for testname, infile, expects in tests:
            if os.path.exists(infile):
                runTest(testname, infile, expects, indexType)

def main():
    # tabixpy.setLogLevel(tabixpy.logging.DEBUG)
    
    t1 = [
        "150.100000",
        "annotated_tomato_150.100000.vcf.gz.tbi",
        # return
        # 1_375_671 - 1_378_902 [-2]
        # 1_392_519 - 1_393_971 [-1]
        # 1_393_980 - 1_395_108
        # 1_393_971
        # 1_395_638 last value
        # 1_395_108
        [
            [0, 1_375_671, None,      1_375_671, 1_395_638, 1_120],
            [0, 1_375_672, None,      1_375_685, 1_395_638, 1_119], #last to last bin
            [0, 1_392_519, None,      1_392_519, 1_395_638,   263], #last bin
            [0, 1_392_520, None,      1_392_520, 1_395_638,   262], #last bin
            [0, 1_395_638, None,      1_395_638, 1_395_638,     1], #last value
            [0, 1_392_520, 1_395_638, 1_392_520, 1_395_632,   261],
            [0, 1_395_639, None,      None,      None,          0],
        ]
    ]

    t2 = [
        "150.SL2.50ch00-01-02",
        "annotated_tomato_150.SL2.50ch00-01-02.vcf.gz",
        [
            [0, 1_375_671,      None, 1_375_671,  1_395_638, 1_120],
            [0, 1_375_672,      None, 1_375_685,  1_395_638, 1_119], #last to last bin
            [0, 1_392_519,      None, 1_392_519,  1_395_638,   263], #last bin
            [0, 1_392_520,      None, 1_392_520,  1_395_638,   262], #last bin
            [0, 1_395_638,      None, 1_395_638,  1_395_638,     1], #last value
            [0, 1_392_520, 1_395_638, 1_392_520,  1_395_632,   261],
            [0, 1_395_639,      None,      None,       None,     0],
            [1,       189,       189,      None,       None,     0],
            [1,       190,       285,       190,        190,     1],
            [1,       190,       286,       190,        190,     1],
            [1,       190,       287,       190,        286,     2]
        ]
    ]

    t3 = [
        "150",
        "annotated_tomato_150.vcf.bgz.tbi",
        [
            [0,  1_375_671,      None,  1_375_671, 21_805_702, 1_107_062],
            [0,  1_375_672,      None,  1_375_685, 21_805_702, 1_107_061], #last to last bin
            [0,  1_392_519,      None,  1_392_519, 21_805_702, 1_106_205], #last bin
            [0,  1_392_520,      None,  1_392_520, 21_805_702, 1_106_204], #last bin
            [0, 21_805_700,      None, 21_805_702, 21_805_702,         1], #last value
            [0,  1_392_520, 1_395_638,  1_392_520,  1_395_632,       261],
            [0, 21_805_703,      None,       None,       None,         0],
            [1,        189,       189,       None,       None,         0],
            [1,        190,       285,        190,        190,         1],
            [1,        190,       286,        190,        190,         1],
            [1,        190,       287,        190,        286,         2]
        ]
    ]
    
    # runTests([
    #     t1,
    #     t2,
    #     t3
    # ], [tabixpy.Formats.TBJ, tabixpy.Formats.TBK])

    # runTests([
    #     t1,
    #     t2,
    #     t3
    # ], [tabixpy.Formats.TBJ])

    # runTests([
    #     t1,
    #     t2,
    # ], [tabixpy.Formats.TBK])

    runTests([
        t1,
        t2,
    ], [tabixpy.Formats.TBJ, tabixpy.Formats.TBK])



    # for p in range(999999999900,1000000000000):
    #     print(p, tabixpy.reg2bin(p, p+1), tabixpy.reg2bins(p, p+1))

    # import gzip
    # print("opening")
    # with open(ingz, "rb") as inf:
    #     for cp, chrom in enumerate(data["names"]):
    #         tabixpy.logger.info("chrom", chrom)

    #         refs  = data["refs"][cp]

    #         n_intv = refs["n_intv"]
    #         tabixpy.logger.info("n_intv", n_intv)
    #         intvs  = refs["intvs"]


            # for intv in intvs:
            #     tabixpy.logger.info(intv)

    #         last_chunk = None
    #         for virtual_offset, chunk_real_begin, chunk_bytes_begin in intvs:
    #             print(virtual_offset, chunk_real_begin, chunk_bytes_begin)
    #             if last_chunk is not None:
    #                 print("    block size", chunk_real_begin - last_chunk)
    #             last_chunk = chunk_real_begin
    #             # inf.seek(0, 0)
    #             inf.seek(chunk_real_begin, 0)
    #             with gzip.open(inf, 'rb') as g:
    #                 b = g.read(tabixpy.BLOCK_SIZE)
    #                 i = b[chunk_bytes_begin:chunk_bytes_begin+20]
    #                 p = int(i.decode().split("\t")[1])
    #                 print("    i", i)
    #                 print("    p", p)


            # n_bin = refs["n_bin"]
            # bins  = refs["bins"]
            # for bin_n, bind in enumerate(bins):
            #     tabixpy.logger.info(f" bin_n {bin_n}")
            #     tabixpy.logger.info(f" chunks {bind['chunks']}")
            #     n_chunk = bind["n_chunk"]
            #     chunk_begins = bind["chunks"]["chunk_begin"]
            #     chunk_ends   = bind["chunks"]["chunk_end"]
            #     for chunk_n in range(n_chunk):
            #         tabixpy.logger.info(f"  chunk_n {chunk_n}")
            #         chunk_begin = chunk_begins[chunk_n]
            #         chunk_end   = chunk_ends[chunk_n]
            #         tabixpy.logger.info(f"   chunk_begin {chunk_begin}")
            #         tabixpy.logger.info(f"   chunk_end   {chunk_end}")
    #                 inf.seek(0, 0)
    #                 inf.seek(chunk_real_begin, 0)
    #                 with gzip.open(inf, 'r') as g:
    #                     b = g.read(tabixpy.BLOCK_SIZE)
    #                     print("    b", b)
    #                     i = b[chunk_bytes_begin:chunk_bytes_begin+20]
    #                     print("    i", i)


        # print("seeking")
        # f.seek(7021611)
        # # f.seek(4631)
        # g=gzip.open(f)
        # d=g.read(16 * 1024)
        # # print("d", d)
        # r=d[4631:4631+100]
        # print(r)
        # # for line in g:
        # #     print("line", line)
        # #     break


if __name__ == "__main__":
    main()