#!/bin/bash

set -xeu

SRC=annotated_tomato_150.vcf.bgz
BN=annotated_tomato_150.SL2.50ch00-01-02.vcf

rm -v ${BN}* || true

tabix -h ${SRC} SL2.50ch00:0-1,395,638 >  ${BN}
tabix    ${SRC} SL2.50ch01:0-2,000,000 >> ${BN}
tabix    ${SRC} SL2.50ch02:0-3,000,000 >> ${BN}

bgzip -c ${BN} > ${BN}.gz

tabix ${BN}.gz

ls -lh ${BN}*
