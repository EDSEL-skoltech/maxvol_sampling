import argparse
import numpy as np
import osgeo.gdal as gdal
import os
from maxvol_cut import rect_maxvol_cut, f_no_cut, f_penal_2D
from tools import norm_data, add_coords, gen_input, extend_score, points_selection, f_no_cut, f_cut_eps, calc_score, good_points_brute_force, idx_to_idx


def data_preparation(wd, data_m, dem_dir):
    fl_names = list(filter(lambda fl: fl.endswith('.tif'), os.listdir(wd+'features/')))
    files = list(map(lambda x: gdal.Open(os.path.join(wd+'features/', x)), fl_names))
    arrays = list(map(lambda x: x.ReadAsArray().flatten(), files))
    shapes = [x.ReadAsArray().shape for x in files]
    nodatas = list(map(lambda x: x.GetRasterBand(1).GetNoDataValue(), files))
    names = list(map(lambda x: x.replace('.tif','').split('.')[0], fl_names))
        
    if dem_dir is None:
        dem_raw = gdal.Open(wd+'/dem.tif')
        dem = dem_raw.ReadAsArray()
        
    else:
        dem_raw = gdal.Open(dem_dir)
        dem = dem_raw.ReadAsArray()
        
    
    dem_flat = dem.flatten()
    dem_nodata = dem_raw.GetRasterBand(1).GetNoDataValue()
    init_dem_shape = dem.shape
    
    #delete nodata
    idx_nodata_0 = np.where(arrays[0] == nodatas[0])[0]
    arrays_no_nodatas = np.zeros((len(arrays[0])-len(idx_nodata_0), len(arrays)))

    idx_dem_nodata = np.where(dem_flat == dem_nodata)[0]
    idx_dem = np.where(arrays[0] != nodatas[0])[0]
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
    
    X = np.vstack((X, X[:1,:]))
    return X, dem_flat, dem_nodata, init_dem_shape, idx_dem

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--max_n_pnts', type=int, default=15, help='number of points', required=True)
    parser.add_argument('--min_n_pnts', type=int, default=15, help='number of points', required=True)
    parser.add_argument('--wd', type=str, help='working directory', required=True)
    parser.add_argument('--dem_dir', default=None, type=str, help='dem directory')
    parser.add_argument('--data_m', type=int, default=0, help='data mode')
    parser.add_argument('--dist_pts', type=float, default=0.3, help='distance between points')
    parser.add_argument('--dist_brd', type=float, default=0.2, help='distance from border')
    
    args = parser.parse_args()
    
    X, dem_flat, dem_nodata, init_dem_shape, idx_dem = data_preparation(args.wd, args.data_m, args.dem_dir)
    X = np.unique(X, axis=1)
    
    
    #function for distance between points
    f_cut = lambda idx, i : f_cut_eps(idx, i, X=X, eps= args.dist_pts)
    #function for distance from border
    f_penal = f_penal_2D(X = X[:, -2], Y = X[:, -1], bnd = args.dist_brd, level = 0.3)

    result = points_selection(X, max_n_pnts = args.max_n_pnts,min_n_pnts = args.min_n_pnts, cut_fun = f_cut, penalty = f_penal) 
    
    #coordinates
    xmin, ymin, xmax, ymax = [37.7928399298814384,51.9023655604198453,37.8064010466814366,51.9077426838598441]

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

    np.savetxt(args.wd+'maxvol_result_'+str(len(result))+'.csv', coord_idx, delimiter=';')
