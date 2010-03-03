#!/bin/bash
echo '<?xml version="1.0"?> <VTKFile type="Collection" version="0.1" ><Collection>' > data.pvd

for file in $(ls *.vtu) ; do 
	echo '<DataSet timestep="0" part="0" file="'${file}'"/>' >> data.pvd
done

#<DataSet timestep="39" part="0" file="temperature000039.vtu"/>



echo '</Collection></VTKFile>' >> data.pvd

