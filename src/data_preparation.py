def data_preparation(wd, data_m, dem_dir):
    """
    Function to orginize tif files in flatten vectos, remove NaN and stack vectors into matrix
    
    """
    fl_names = list(filter(lambda fl: fl.endswith('.tif'), os.listdir(wd+'/local_tifs/')))
    files = list(map(lambda x: gdal.Open(os.path.join(wd+'/local_tifs/', x)), fl_names))
    #     files = list(map(lambda x: gdal.Open(os.path.join(wd, x)), fl_names))
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



    idx_nodata_0 = np.where(dem_flat == dem_nodata)[0]


    arrays_no_nodatas = np.zeros((len(arrays[0])-len(idx_nodata_0), len(arrays)))

    idx_dem_nodata = np.where(dem_flat == dem_nodata)[0]
    idx_dem = np.where(dem_flat != dem_nodata)[0]
    print(idx_dem.shape)
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
    return X, dem_flat, dem_nodata, init_dem_shape, idx_dem
