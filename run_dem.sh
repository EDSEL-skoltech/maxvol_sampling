#!/bin/bash
for i in {1..12} 
do
	bash run_maxvol.sh --path_to_file ./tiff_urupinsk/field_$i.tif --path_to_output ./tiff_urupinsk/result_features_$i --number_of_points 19
done
