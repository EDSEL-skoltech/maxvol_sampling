#!/bin/bash 
echo There were $# DEM files!
echo Using meathod \$@: $@
if [[ $# -eq "0" ]]; then 
	echo "Hey, man! I need TIF DEM to compute MaxVol, yo!";
	else

		for param in "$@"
		do 
			for j in 25 30
			do
				for i in $(seq 0 1000)
				do
					echo $i
					cp $param dem.tif
					rm -f ./features/*tif
					python noisy_dem.py --path_to_dem ./dem.tif --path_to_noisy_result ./noisy_dem.tif
					docker run --rm -it -v $PWD:/data whitebox_tools --run=Aspect -i=noisy_dem.tif -o=aspect.tif
					docker run --rm -it -v $PWD:/data whitebox_tools --run=Slope -i=noisy_dem.tif -o=slope.tif
					docker run --rm -it -v $PWD:/data whitebox_tools --run=WetnessIndex --sca=noisy_dem.tif --slope='slope.tif' -o=wetnessindex.tif
					mv aspect.tif slope.tif wetnessindex.tif ./features 
					cp ../../../DEM_files/ndvi_features/sink.tif ./features
					cp ../../../DEM_files/ndvi_features/ndvi_dem.tif ./features
					python maxvol_points_selection_with_whitebox_noisy.py --max_n_pnts $j --min_n_pnts $j --wd $PWD --dist_pts 0.3
					python clhs_noisy.py  --wd $PWD --n_pnts $j 
					rm -f ./noisy_dem.tif dem.tif
				done
			done
		done	
fi
