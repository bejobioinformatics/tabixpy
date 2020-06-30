set -xeu

rm -v *.tbj *.log || true
rm -v __pycache__ || true

pypy test.py  2>&1 > test.log
