#!/bin/bash

j=1

for i in $(cat problems.list) ; do 
	DIR=pol${j}_${i}	
	mkdir $DIR
	cd $DIR
	cmake -DPOLORDER:STRING="${1}" -DGRIDTYPE:STRING="${2}" -DPROBLEM:STRING="${i}"
		-DCMAKE_VERBOSE_MAKEFILE:STRING="OFF" \
		-DCXX_FLAGS:STRING="-DNDEBUG -O3 -fomit-frame-pointer -funroll-loops -w" ../../dune-stokes/
	ln -s ../grid_2d.dgf
	ln -s ../gdbinit .gdbinit
	ln -s ../test.param 
	mkdir data
	cd data
	ln -s ../../eoc-template.tex
	cd ..
done
