#!/bin/bash

cat header_template.tex > eoc.tex

for i in $(ls); do
	if [ -d "${i}" ] ; then
		echo "\\input{${i}/data/eoc-file_body}" >> eoc.tex
	fi
done

cat footer_template.tex >> eoc.tex

pdflatex eoc.tex
