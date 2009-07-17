#!/bin/bash

for i in $(ls) ; do
        if [ -d "$i" ] ; then
                cd $i
                screen -S "${i}" -d -m ./dune_stokes test.param
                cd ..
        fi
done

