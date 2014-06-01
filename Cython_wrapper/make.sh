#!/bin/bash
rm -rf build FMIndex.cpp FMIndex.so
export CC="g++ -std=c++11"
export ARCHFLAGS="-arch x86_64"
python setup.py build_ext --inplace
