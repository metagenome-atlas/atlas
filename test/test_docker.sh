#! /bin/bash

pip install -e .

atlas --help

cd test

./test_local.sh #--conda-frontend conda