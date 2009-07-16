#!/bin/bash

if [ $3 == 1 ] ; then
	DIR=pol${1}_bfg_${2}
else
	DIR=pol${1}_no-bfg_${2}
fi

mkdir $DIR
cd $DIR
cmake -DPOLORDER:STRING="${1}" -DGRIDTYPE:STRING="${2}" -DCMAKE_VERBOSE_MAKEFILE:STRING="OFF" \
	-DCXX_FLAGS:STRING="-DNDEBUG -O3 -fomit-frame-pointer -funroll-loops -w" ../../dune-stokes/
ln -s ../grid_2d.dgf
ln -s ../gdbinit .gdbinit
cp ../test.param .
echo "do-bfg: $3" >> test.param
mkdir data
ln -s ../eoc-template.tex data/eoc-template.tex
make

