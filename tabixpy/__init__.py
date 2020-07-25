import os
from .tabixpy import *

absp = os.path.abspath(os.path.dirname(__file__))
print("absp", absp)

__author__       = "Saulo Aflitos"
__author_email__ = "saulobejo@users.noreply.github.com"
__url__          = "https://github.com/bejobioinformatics/tabixpy"
__description__  = "Tabix reader written 100% in Python"
__package_name__ = "tabixpy"
___keywords__    = "Tabix Genomics VCF SNP"

with open(os.path.join(absp, "README.md"), "rt") as fh:
    __long_description__ = fh.read()

with open(os.path.join(absp, "LICENSE.txt"), "rt") as fh:
    __license__ = fh.read()

with open(os.path.join(absp, "VERSION"), "rt") as fh:
    __version__ = fh.read()
