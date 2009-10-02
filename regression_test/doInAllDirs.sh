#!/bin/bash

for i in $(find . -mindepth 1  -maxdepth 1 -type d ) ; do 
	cd $i 
	${@}  
	cd ..  
done

