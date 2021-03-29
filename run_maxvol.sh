#!/bin/bash

CLEAR='\033[0m'
RED='\033[0;31m'

function usage() {
  if [ -n "$1" ]; then
    echo -e "${RED}👉 $1${CLEAR}\n";
  fi
  echo "Usage: $0 [-f path_to_file] [-d path_to_output]"
  echo "  -f, --path_to_file   Path to DEM file"
  echo "  -d, --path_to_output  Path to new folder with terrain features"
  echo "  -n, --number_of_points  Number of sampling points"
  echo "  -v, --verbose  Verbose"
  echo "  -h, --help  Help"
  echo ""
  echo "Example: $0 --path_to_file ./DEM.tif --path_to_output ./terrain --number_of_points 10"
  exit 1
}

# parse params
while [[ "$#" > 0 ]]; do case $1 in
  -f|--path_to_file) PATH_TO_DEM="$2"; shift;shift;;
  -d|--path_to_output) PATH_TO_DIR="$2";shift;shift;;
  -n|--number_of_points) NUMBER_OF_POINTS="$2";shift;shift;;
  -v|--verbose) VERBOSE=1;shift;;
  -h|--help) usage;shift;;
  *) usage "Unknown parameter passed: $1"; shift; shift;;
esac; done

# verify params
if [ -z "$PATH_TO_DEM" ]; then usage "Input DEM file path is not set"; fi;
if [ -z "$PATH_TO_DIR" ]; then usage "Path to output is not set."; fi;
if [ -z "$NUMBER_OF_POINTS" ]; then usage "Number of points is not set."; fi;

# Test folder and files in folder; if OK -> mkdir 
if [ -d $PATH_TO_DIR ]; then
    if [ ! -z "$(ls -A $PATH_TO_DIR)" ]; then
        echo "Dir is not empty"
        exit
    else
        echo "Output path - Ok"
    fi
else
    mkdir $PATH_TO_DIR
fi

cp $PATH_TO_DEM ./wb_dem.tif

docker run --rm -it -v $PWD:/data mishagrol/whitebox_tools --run=Aspect -i=wb_dem.tif -o=aspect.tif
docker run --rm -it -v $PWD:/data mishagrol/whitebox_tools --run=Slope -i=wb_dem.tif -o=slope.tif
docker run --rm -it -v $PWD:/data mishagrol/whitebox_tools --run=WetnessIndex --sca=wb_dem.tif --slope='slope.tif' -o=wetnessindex.tif
docker run --rm -it -v $PWD:/data mishagrol/whitebox_tools --run=Sink -i=wb_dem.tif -o=sink.tif
mv aspect.tif slope.tif wetnessindex.tif sink.tif ./$PATH_TO_DIR

python3 maxvol_sampling.py --max_n_pnts $NUMBER_OF_POINTS --min_n_pnts $NUMBER_OF_POINTS --path_to_DEM ./wb_dem.tif --path_to_DEM_features ./$PATH_TO_DIR --dist_pts 0.1
