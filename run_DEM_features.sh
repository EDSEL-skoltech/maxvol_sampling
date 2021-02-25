#!/bin/bash

# usage function
function usage()
{
   echo "
   Usage: $0 [--path_to_file DEM.tif] [--path_to_dir ./dem_features] [--verbose] [--help]
   optional arguments:
     -h, --help           show this help message and exit
     -f, --path_to_file   path to DEM file
     -p, --path_to_dir    path to output folder with terrain features
     -v, --verbose        increase the verbosity of the bash script"
}  

# initialize variables
progname=$(basename $0)
verbose=0
path_to_dir="./dem_features"
path_to_file="./demdegrdr.tif"
number_of_points=10
time_str=0

# use getopt and store the output into $OPTS
# note the use of -o for the short options, --long for the long name options
# and a : for any option that takes a parameter
OPTS=$(getopt -o "hf:p:n:v" --long "help,path_to_file:,path_to_dir:,number_of_points:,verbose," -n "$progname" -- "$@")
if [ $? != 0 ] ; then echo "Error in command line arguments." >&2 ; usage; exit 1 ; fi
eval set -- "$OPTS"

while true; do
  # uncomment the next line to see how shift is working
  echo "\$1:\"$1\" \$2:\"$2\""
  case "$1" in
    -h | --help ) usage; exit; ;;
    -f | --path_to_file ) path_to_file="$2"; shift 2 ;;
    -p | --path_to_dir ) path_to_dir="$2"; shift 2 ;;
    # -n | --number_of_points) number_of_points="$2"; shift 1 ;;
    -v | --verbose ) verbose=$((verbose + 1)); shift ;;
    -- ) shift; break ;;
    * ) break ;;
  esac
done

if (( $verbose > 0 )); then

   # print out all the parameters we read in
   cat <<EOM
   path_to_file=$path_to_file
   path_to_dir=$path_to_dir
   number_of_points=$number_of_points
   verbose=$verbose
EOM
fi

# The rest of your script below
if [ -d $path_to_dir ]; then
    if [ ! -z "$(ls -A $path_to_dir)" ]; then
        echo "Dir is not empty"
        exit
    else
        echo "Ok"
    fi
else
    mkdir $path_to_dir
fi
# echo $path_to_file
# echo $path_to_dir

cp $path_to_file ./wb_dem.tif

echo $number_of_points

# docker run --rm -it -v $PWD:/data mishagrol/whitebox_tools --run=Aspect -i=wb_dem.tif -o=aspect.tif
# docker run --rm -it -v $PWD:/data mishagrol/whitebox_tools --run=Slope -i=wb_dem.tif -o=slope.tif
# docker run --rm -it -v $PWD:/data mishagrol/whitebox_tools --run=WetnessIndex --sca=wb_dem.tif --slope='slope.tif' -o=wetnessindex.tif
# docker run --rm -it -v $PWD:/data mishagrol/whitebox_tools --run=Sink -i=wb_dem.tif -o=sink.tif

# mv aspect.tif slope.tif wetnessindex.tif sink.tif -t ./$path_to_dir

# python3 maxvol_sampling.py --max_n_pnts $number_of_points --min_n_pnts $number_of_points --path_to_DEM ./wb_dem.tif --path_to_DEM_features ./$path_to_dir --dist_pts 0.1

