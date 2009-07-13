#!/bin/bash

mkdir build-pol${1}
cd build-pol${1}
cmake -DPOLORDER="${1}" ..
cp -r ../src/{grid_2d.dgf,test.param} .
mkdir data
cp ../src/data/eoc* data/
 
