import sys

sys.path.insert(0, '..')

import tabixpy

if __name__ == "__main__":
    infile     = sys.argv[1]
    ingz, inid = tabixpy.get_filenames(infile)

    # tabixpy.logger.setLevel(tabixpy.logging.DEBUG)

    data       = tabixpy.read_tabix(ingz)

    tabixpy.save(data, ingz, compress=True)
