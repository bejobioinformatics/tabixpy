set -xeu

rm -v *.tbj *.log || true
rm -v __pycache__ || true

pypy test.py annotated_tomato_150.100000.vcf.gz.tbi 2>&1 > annotated_tomato_150.100000.vcf.gz.tbj.log

# pypy test.py annotated_tomato_150.vcf.bgz.tbi 2>&1 > annotated_tomato_150.vcf.bgz.tbj.log
