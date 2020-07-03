#!/bin/bash

set -xeu

rm -v *.tbj *.log || true
rm -v *.tbj.json  || true
rm -v __pycache__ || true

python3 test.py  2>&1 > test.log

(cd ..; python -c 'import tabixpy; tabixpy.genAllBlocks("tests/annotated_tomato_150.100000.vcf.gz")')
(cd ..; python -c 'import tabixpy; _= tabixpy.loadAllBlocks("tests/annotated_tomato_150.100000.vcf.gz")')

(cd ..; python -c 'import tabixpy; tabixpy.genAllBlocks("tests/annotated_tomato_150.SL2.50ch00-01-02.vcf.gz")')
(cd ..; python -c 'import tabixpy; _= tabixpy.loadAllBlocks("tests/annotated_tomato_150.SL2.50ch00-01-02.vcf.gz")')


# f=annotated_tomato_150.100000.vcf.gz.tbj
# f=annotated_tomato_150.SL2.50ch00-01-02.vcf.gz.tbj
f=annotated_tomato_150.vcf.bgz.tbj
gunzip -kc $f > $f.json
