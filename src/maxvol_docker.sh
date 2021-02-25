#!/bin/bash 
echo There were $# DEM files!
echo Using meathod \$@: $@
conda activate maxvol_soil_sampling
if [[ $# -eq "0" ]]; then 
	echo "Hey, man! I need TIF DEM to compute MaxVol, yo!";
	else

		for param in "$@"
		do 
			cp $param dem.tif
			rm -f ./features/*tif
			docker run --rm -it -v $PWD:/data whitebox-tools --run=Aspect -i=$param -o=aspect.tif
			docker run --rm -it -v $PWD:/data whitebox-tools --run=Slope -i=$param -o=slope.tif
			docker run --rm -it -v $PWD:/data whitebox-tools --run=WetnessIndex --sca=$param --slope='slope.tif' -o=wetnessindex.tif
			docker run --rm -it -v $PWD:/data whitebox-tools --run=Sink -i=$param -o=sink.tif
			mv aspect.tif slope.tif wetnessindex.tif sink.tif ./features  
			python maxvol_points_selection_with_whitebox_grol.py --max_n_pnts 10 --min_n_pnts 10 --wd $PWD --dist_pts 0.25
		done	
fi
