#!/bin/bash

NUMBER_OF_POINTS=15
cp dem.tif wb_dem.tif

PATH_TO_DIR=./urupinsk/

docker run --rm -it -v $PWD:/data mishagrol/whitebox_tools --run=Aspect -i=wb_dem.tif -o=aspect.tif
docker run --rm -it -v $PWD:/data mishagrol/whitebox_tools --run=Slope -i=wb_dem.tif -o=slope.tif
docker run --rm -it -v $PWD:/data mishagrol/whitebox_tools --run=WetnessIndex --sca=wb_dem.tif --slope='slope.tif' -o=wetnessindex.tif
docker run --rm -it -v $PWD:/data mishagrol/whitebox_tools --run=Sink -i=wb_dem.tif -o=sink.tif

mv aspect.tif slope.tif wetnessindex.tif sink.tif ./$PATH_TO_DIR

python3 maxvol_sampling.py --max_n_pnts $NUMBER_OF_POINTS --min_n_pnts $NUMBER_OF_POINTS --path_to_DEM ./wb_dem.tif --path_to_DEM_features ./$PATH_TO_DIR --dist_pts 0.1
