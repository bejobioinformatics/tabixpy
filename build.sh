#!/bin/bash

#https://packaging.python.org/tutorials/packaging-projects/

set -xeu

rm -rfv build/ dist/ tabixpy.egg-info/ tabixpy-* || true

#python3 -m pip install --user --upgrade setuptools wheel
#python3 -m pip install --user --upgrade twine

python3 setup.py sdist bdist_wheel

twine check dist/tabixpy-*-py3-none-any.whl
