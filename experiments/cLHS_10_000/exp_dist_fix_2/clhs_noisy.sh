#!/bin/bash 
echo There were $# DEM files!
echo Using meathod \$@: $@
if [[ $# -eq "0" ]]; then 
	echo "Hey, man! I need TIF DEM to compute MaxVol, yo!";
	else

		for param in "$@"
		do 
			for i in $(seq 0 500)
			do
				echo $i
				cp $param dem.tif
				rm -f ./features/*tif
				rm -f ./noisy_dem.tif
				python noisy_dem.py --path_to_dem ./dem.tif --path_to_noisy_result ./noisy_dem.tif
				docker run --rm -it -v $PWD:/data whitebox-tools --run=Aspect -i=noisy_dem.tif -o=aspect.tif
				docker run --rm -it -v $PWD:/data whitebox-tools --run=Slope -i=noisy_dem.tif -o=slope.tif
				docker run --rm -it -v $PWD:/data whitebox-tools --run=WetnessIndex --sca=noisy_dem.tif --slope='slope.tif' -o=wetnessindex.tif
				#docker run --rm -it -v $PWD:/data whitebox-tools --run=Sink -i=noisy_dem.tif -o=sink.tif
				mv aspect.tif slope.tif wetnessindex.tif ./features  
				cp ./sink_file/sink.tif ./features
				python clhs_noisy.py  --wd $PWD --n_pnts 15 
				rm dem.tif
			done
		done	
fi
