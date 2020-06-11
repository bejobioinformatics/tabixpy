set -xeu

rm -v *.tbj || true
rm -v __pycache__ || true

pypy test.py annotated_tomato_150.100000.vcf.gz.tbi

# pypy test.py annotated_tomato_150.vcf.bgz.tbi
