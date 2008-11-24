#!/bin/bash
if [ $# -lt 1 ] ; then
   echo 'setting up dune'
   ./dune-common/bin/dunecontrol --opts=config.opts.wwu_no_documentation all
   exit 0 
fi
echo 'setting up ' $@
./dune-common/bin/dunecontrol --only=$@ --opts=config.opts.wwu_no_documentation all
exit 0

