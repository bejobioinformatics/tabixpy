import sys

sys.path.insert(0, '..')

import tabixpy

if __name__ == "__main__":
    infile     = sys.argv[1]
    ingz, inid = tabixpy.get_filenames(infile)

    tabixpy.logger.setLevel(tabixpy.logging.DEBUG)

    data       = tabixpy.read_tabix(ingz)

    tabixpy.save(data, ingz, compress=False)


    # for p in range(999999999900,1000000000000):
    #     print(p, tabixpy.reg2bin(p, p+1), tabixpy.reg2bins(p, p+1))

    # import gzip
    # print("opening")
    # with open(ingz, "rb") as inf:
    #     for cp, chrom in enumerate(data["names"]):
    #         print("chrom", chrom)

    #         refs  = data["refs"][cp]

    #         n_intv = refs["n_intv"]
    #         print("n_intv", n_intv)
    #         intvs  = refs["intvs"]

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


    #         n_bin = refs["n_bin"]
    #         bins  = refs["bins"]
    #         for bin_n, bind in enumerate(bins):
    #             print(" bin_n", bin_n)
    #             n_chunk = bind["n_chunk"]
    #             chunk_real_begins  = bind["chunks"]["chunk_real_begin"]
    #             chunk_bytes_begins = bind["chunks"]["chunk_bytes_begin"]
    #             for chunk_n in range(n_chunk):
    #                 print("  chunk_n", chunk_n)
    #                 chunk_real_begin  = chunk_real_begins[chunk_n]
    #                 chunk_bytes_begin = chunk_bytes_begins[chunk_n]
    #                 print("   chunk_real_begin ", chunk_real_begin)
    #                 print("   chunk_bytes_begin", chunk_bytes_begin)
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