#!/bin/bash
if [ $# -lt 1 ] ; then
   echo 'setting up dune'
   ./dune-common/bin/dunecontrol --opts=config.opts
   exit 0 
fi
echo 'setting up ' $@
./dune-common/bin/dunecontrol --only=$2 --opts=config.opts all
exit 0

