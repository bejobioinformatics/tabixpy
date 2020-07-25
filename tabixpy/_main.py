import os
import sys
import argparse

#sys.path.insert(0, '..')

import tabixpy

class PrinterAction(argparse.Action):
    # VAR = "NONE"

    def __init__(self, option_strings, dest, nargs=None, **kwargs):
        if nargs is not None:
            raise ValueError("nargs not allowed")
        super(PrinterAction, self).__init__(option_strings, dest, nargs=0, **kwargs)
    def __call__(self, parser, namespace, values, option_string=None):
        print(self.VAR)
        sys.exit(0)
        #print('%r %r %r' % (namespace, values, option_string))
        #setattr(namespace, self.dest, values)

class DescriptionVar:
    __slot__ = True
    VAR = tabixpy.__long_description__

class LicenseVar:
    __slot__ = True
    VAR = tabixpy.__license__

class Description(PrinterAction, DescriptionVar):
    pass

class License(PrinterAction, LicenseVar):
    pass

parser = argparse.ArgumentParser(
    description=f'Read Tabix index file and create either TBJ or TBK indexes. {sys.argv[0]} Version {tabixpy.__version__}',
    prog=sys.argv[0],
    usage='%(prog)s [options] <FORMAT> <INPUT FILE>'
    )

parser.add_argument('--version', action='version', version=f'%(prog)s {tabixpy.__version__}', help='Print version and exit')

parser.add_argument('--description', '-d', action=Description, help='Print description and exit')

parser.add_argument('--license', '-l', action=License, help='Print license and exit')

parser.add_argument(
    '--no-overwrite'     if tabixpy.TabixDefaults.OVERWRITE else '--overwrite',
    '-noo'               if tabixpy.TabixDefaults.OVERWRITE else '-o',
    dest="overwrite",
    action='store_false' if tabixpy.TabixDefaults.OVERWRITE else 'store_true',
    help  =f'Overwrite file if exists. default: {tabixpy.TabixDefaults.OVERWRITE}'
)

parser.add_argument(
    '--no-compress'      if tabixpy.TabixDefaults.COMPRESS else '--compress',
    '-noc'               if tabixpy.TabixDefaults.COMPRESS else '-c',
    dest="compress",
    action='store_false' if tabixpy.TabixDefaults.COMPRESS else 'store_true',
    help  =f'Compress TBJ file. default: {tabixpy.TabixDefaults.COMPRESS}'
)

parser.add_argument('--verbose', '-v', action='count', default=0, help='Verbosity level. accepts multiple.')

fmts = list(tabixpy.Formats.__members__.keys())
parser.add_argument('format', metavar='FORMAT'    , choices=fmts, default=tabixpy.Formats.DEFAULT, help=f'Index format. options: {", ".join(fmts)}. default: {tabixpy.Formats.DEFAULT.name}')

parser.add_argument('infile', metavar='INPUT FILE', nargs=1, help='VCF file to index: *.vcf.[b]gz[.tbi] file')

def main():
    parsed    = parser.parse_args()

    indexType   = tabixpy.Formats[parsed.format]
    verbosity   = parsed.verbose * 10
    infile      = parsed.infile[0]
    compress    = parsed.compress
    overwrite   = parsed.overwrite
    description = parsed.description
    license     = parsed.license

    if not os.path.exists(infile):
        print(f"Input file '{infile}' does not exists.\n-----------------", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    print(f"""
compress  = {compress}
overwrite = {overwrite}
verbosity = {verbosity}
indexType = {indexType}
infile    = {infile}
""")

    tabix  = tabixpy.Tabix(infile, indexType=indexType, logLevel=verbosity)
    tabix.create(overwrite=overwrite, compress=compress)

if __name__ == "__main__":
    main()
