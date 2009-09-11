#!/bin/bash

for i in $(find . -mindepth 1  -maxdepth 1 -type d ) ; do
        cd $i
        make && ./dune_stokes test.param || jabnotif "compile failed on $HOSTNAME, $PWD" rene@graasmilk.net
        cd ..
done

jabnotif "regression test done $HOSTNAME, $PWD" rene@graasmilk.net
























