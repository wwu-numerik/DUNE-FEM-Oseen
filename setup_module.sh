#!/bin/bash

if [ x$2 = x ] ; then
	./dune-common/bin/dunecontrol --module=$1 --opts=config.opts.wwu_no_documentation all
else
	./dune-common/bin/dunecontrol --module=$1 --opts=${2} all	
fi
