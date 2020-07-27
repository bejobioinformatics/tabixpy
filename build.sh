#!/bin/bash

#https://packaging.python.org/tutorials/packaging-projects/

set -xeu

if "1" == "1"; then
    cd tabixpy/tests

    bash run.sh

    rm -v *.tbj *.tbk *.log *.json

    cd -
fi


rm -rfv build/ dist/ tabixpy.egg-info/ tabixpy-* || true

#python3 -m pip install --user --upgrade setuptools wheel
#python3 -m pip install --user --upgrade twine

python3 setup.py sdist bdist_wheel

twine check dist/tabixpy-*-py3-none-any.whl
