#!/bin/bash
cd dune-common
make clean
make distclean
cd ../dune-grid
make clean
make distclean
cd ../dune-istl
make clean
make distclean
cd ../dune-fem
make clean
make distclean
cd ../dune-grid-howto
make clean
make distclean
cd ../dune-fem-howto
make clean
make distclean
cd ../mirko-stokes
make clean
make distclean
cd ../dune-stokes
make clean
make distclean


