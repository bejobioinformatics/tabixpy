#!/bin/bash

set -xeu
set -o pipefail

rm -v __pycache__ || true
rm -v *.log       || true
rm -v *.tbj       || true
rm -v *.tbj.json  || true
rm -v *.tbk       || true

PROG=python3
#PROG=pypy3

${PROG} test.py 2>&1 > test.log

(cd ..; time ${PROG} -c 'import tabixpy; _= tabixpy.genVcfGzPy("tests/annotated_tomato_150.100000.vcf.gz")') | tee test.log
(cd ..; time ${PROG} -c 'import tabixpy; _= tabixpy.loadVcfGzPy("tests/annotated_tomato_150.100000.vcf.gz")') | tee test.log

(cd ..; time ${PROG} -c 'import tabixpy; _= tabixpy.genVcfGzPy("tests/annotated_tomato_150.SL2.50ch00-01-02.vcf.gz")') | tee test.log
(cd ..; time ${PROG} -c 'import tabixpy; _= tabixpy.loadVcfGzPy("tests/annotated_tomato_150.SL2.50ch00-01-02.vcf.gz")') | tee test.log

(cd ..; time ${PROG} -c 'import tabixpy; _= tabixpy.genVcfGzPy("tests/annotated_tomato_150.vcf.bgz")') | tee test.log
(cd ..; time ${PROG} -c 'import tabixpy; _= tabixpy.loadVcfGzPy("tests/annotated_tomato_150.vcf.bgz")') | tee test.log


# f=annotated_tomato_150.100000.vcf.gz.tbj
# f=annotated_tomato_150.SL2.50ch00-01-02.vcf.gz.tbj
if [[ -f "annotated_tomato_150.vcf.bgz.tbj" ]]; then
    f=annotated_tomato_150.vcf.bgz.tbj
    gunzip -kc $f > $f.json
fi
