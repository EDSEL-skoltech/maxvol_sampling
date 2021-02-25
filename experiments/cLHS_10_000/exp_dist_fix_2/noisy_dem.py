import numpy as np
import argparse
import os
from osgeo import gdal
from osgeo.gdalconst import GDT_Float32
import sys


def noisy_dem(raster_input, raster_output):
    in_data, out_data = None, None
    # raster_input = './src/10_3857.tif'
    # open input raster
    in_data = gdal.Open(raster_input)
    if in_data is None:
        print ('Unable to open %s' % raster_input)
        #return None

    # read in data from first band of input raster
    band1 = in_data.GetRasterBand(1)
    rows = in_data.RasterYSize
    cols = in_data.RasterXSize
    vals = band1.ReadAsArray(0, 0, cols, rows)

    # raster_output = './src/test_noisy.tif'

    driver = in_data.GetDriver()
    out_data = driver.Create(raster_output, cols, rows, 1, GDT_Float32)

    dem_data = np.array(vals)

    #dem_data[dem_data < 0] = -32767.
    std_of_DEM = np.std(dem_data[dem_data>0])
    reduced_std = 0.01*std_of_DEM
    # reduced_std = 0.
    print('std for noisy:', reduced_std)

    mask = np.where(dem_data>0)
    dem_data[mask] += np.random.normal(0, reduced_std, size=dem_data[mask].shape) # plus norm dist
    # dem_data[mask] += 10. # plus arbitrary value
    dem_data[dem_data < 0] = -32767.
    out_band = out_data.GetRasterBand(1)
    out_band.WriteArray(dem_data)
    out_band.FlushCache()
    out_band.SetNoDataValue(-32767.)

    out_data.SetGeoTransform(in_data.GetGeoTransform())
    out_data.SetProjection(in_data.GetProjection())
    del out_data
    return raster_output

if __name__=='__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--path_to_dem', type=str, help='path to DEM.tif file', required=True)
    parser.add_argument('--path_to_noisy_result', type=str, help='path to generated noisy DEM file', required=True)
    args = parser.parse_args()

    noisy_dem(args.path_to_dem, args.path_to_noisy_result)
