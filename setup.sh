#!/bin/bash

if [ x$2 = x ] ; then
	./dune-common/bin/dunecontrol --opts=config.opts.wwu_no_documentation all
else
	./dune-common/bin/dunecontrol --opts=${2} all	
fi
