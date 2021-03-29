import argparse
import numpy as np
import osgeo.gdal as gdal
import os
# from maxvol_cut import rect_maxvol_cut, f_no_cut, f_penal_2D
# from tools import norm_data, add_coords, gen_input, extend_score, points_selection, f_no_cut, f_cut_eps, calc_score, good_points_brute_force, idx_to_idx
import subprocess
import csv

from src.util import add_coords, gen_input, points_selection, f_no_cut, f_cut_eps, idx_to_idx
from src.util import rect_maxvol_cut, f_no_cut, f_penal_2D
def data_preparation(path_to_dem, path_to_features, data_m, dem_dir):
    """
    Function to orginize tif files in flatten vectos, remove NaN and stack vectors into matrix
    
    """
    fl_names = list(filter(lambda fl: fl.endswith('.tif'), os.listdir(path_to_features)))
    files = list(map(lambda x: gdal.Open(os.path.join(path_to_features, x)), fl_names))
    #     files = list(map(lambda x: gdal.Open(os.path.join(wd, x)), fl_names))
    arrays = list(map(lambda x: x.ReadAsArray().flatten(), files))
    shapes = [x.ReadAsArray().shape for x in files]
    nodatas = list(map(lambda x: x.GetRasterBand(1).GetNoDataValue(), files))
    names = list(map(lambda x: x.replace('.tif','').split('.')[0], fl_names))

    if dem_dir is None:
        dem_raw = gdal.Open(path_to_dem)
        dem = dem_raw.ReadAsArray()

    else:
        dem_raw = gdal.Open(dem_dir)
        dem = dem_raw.ReadAsArray()

    xmin, xres, xskew, ymax, yskew, yres  = dem_raw.GetGeoTransform()
    xmax = xmin + (dem_raw.RasterXSize * xres)
    ymin = ymax + (dem_raw.RasterYSize * yres)

    boundary_box = {'xmin':xmin, 'xmax':xmax, 'ymin':ymin, 'ymax':ymax}

    dem_flat = dem.flatten()
    dem_nodata = dem_raw.GetRasterBand(1).GetNoDataValue()
    init_dem_shape = dem.shape



    idx_nodata_0 = np.where(dem_flat == dem_nodata)[0]


    arrays_no_nodatas = np.zeros((len(arrays[0])-len(idx_nodata_0), len(arrays)))

    idx_dem_nodata = np.where(dem_flat == dem_nodata)[0]
    idx_dem = np.where(dem_flat != dem_nodata)[0]
    # print(idx_dem.shape)
    dem_no_nodata = np.delete(dem_flat, idx_dem_nodata)


    for i in range(len(arrays)):
        idx_nodata = np.where(arrays[i] == nodatas[i])[0]
        array = arrays[i].copy()
        array[idx_nodata]=0

        arrays_no_nodatas[:,i]  = np.delete(array, idx_nodata_0)
    data_arr = arrays_no_nodatas.copy()

    # Prepare data
    # U can normilize data, and/or add coords to it
    mode = data_m # Change to 0, 1, 2 or 3
    X, fn_X_embedded = gen_input(mode, data_arr, shapes, idx_dem)
#     X = np.vstack((X, X[:1,:]))
    return X, dem_flat, dem_nodata, init_dem_shape, idx_dem, boundary_box

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--max_n_pnts', type=int, default=15, help='number of points', required=True)
    parser.add_argument('--min_n_pnts', type=int, default=15, help='number of points', required=True)
    # parser.add_argument('--wd', type=str, help='working directory', required=True)
    parser.add_argument('--path_to_DEM_features', type=str, help='path to folder with terrain features', required=True)
    parser.add_argument('--path_to_DEM', type=str, help='path to DEM', required=True)
    parser.add_argument('--dem_dir', default=None, type=str, help='dem directory')
    parser.add_argument('--data_m', type=int, default=3, help='data mode')
    parser.add_argument('--dist_pts', type=float, default=0.1, help='distance between points')
    parser.add_argument('--dist_brd', type=float, default=0.2, help='distance from border')
    
    args = parser.parse_args()

    
    X, dem_flat, dem_nodata, init_dem_shape, idx_dem, boundary_box = data_preparation(args.path_to_DEM, args.path_to_DEM_features, args.data_m, args.dem_dir)
    
    print('Selecting points...')
    #function for distance between points
    f_cut = lambda idx, i : f_cut_eps(idx, i, X=X, eps= args.dist_pts)
    #function for distance from border
    f_penal = f_penal_2D(X = X[:, -2], Y = X[:, -1], bnd = args.dist_brd, level = 0.3)

    result = points_selection(X, max_n_pnts = args.max_n_pnts,min_n_pnts = args.min_n_pnts, cut_fun = f_cut, penalty = f_penal) 
    
    #coordinates
    # xmin, ymin, xmax, ymax = [37.7928399,51.90236556, 37.8064010,51.90774268]
    
    xmin = boundary_box['xmin']
    xmax = boundary_box['xmax']
    ymin = boundary_box['ymin']
    ymax = boundary_box['ymax']


    dem_flat_img = dem_flat.copy()-np.min(dem_flat)
    dem_flat_img[np.where(dem_flat == dem_nodata)] = float('NaN')
    st = dem_flat_img.reshape(init_dem_shape)

    xv = np.linspace(xmin,xmax, num=st.shape[1])
    yv = np.linspace(ymax,ymin, num=st.shape[0])
    coords = np.meshgrid(xv,yv)

    mask = idx_dem


    #select corresponding points by indecies
    y_c,x_c = coords[0].flatten()[mask, None],coords[1].flatten()[mask, None]
    y_idx, x_idx = y_c[result],x_c[result]
    coord_idx = np.hstack((y_idx,x_idx))
    print('Done!')
    np.savetxt('maxvol_result_'+str(len(result))+'.csv', coord_idx, delimiter=';')
