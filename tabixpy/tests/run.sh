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

for f in *.tbj; do
    if [[ -f "${f}" ]]; then
        gunzip -kc $f > $f.json
    fi
done
