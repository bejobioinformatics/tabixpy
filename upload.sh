#!/bin/bash

set -xeu

twine check dist/tabixpy-*-py3-none-any.whl

# python3 -m twine upload --repository pypi dist/*

python3 -m twine upload --repository tabixpy dist/*
