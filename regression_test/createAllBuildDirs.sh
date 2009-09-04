#!/bin/bash

for pol in $(cat polorder.list) ; do 
	for grid in $(cat grids.list) ; do 
		for prob in $(cat problems.list) ; do 
			DIR=pol${pol}_${prob}_${grid}
			mkdir $DIR
			cd $DIR
			cmake -DPOLORDER:STRING="${pol}" -DGRIDTYPE:STRING="${grid}" -DPROBLEM:STRING="${prob}" \
				-DCMAKE_VERBOSE_MAKEFILE:STRING="OFF" \
				-DCXX_FLAGS:STRING="-DNDEBUG -O3 -fomit-frame-pointer -funroll-loops -w" ../../dune-stokes/
			ln -s ../grid_2d.dgf
			ln -s ../gdbinit .gdbinit
			ln -s ../test.param 
			mkdir data
			cd data
			ln -s ../../eoc-template.tex
			cd ../..
		done
	done
done
