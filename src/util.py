
import numpy as np
import argparse
import osgeo.gdal as gdal
from scipy.spatial import voronoi_plot_2d, Voronoi
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull, convex_hull_plot_2d
from numpy import genfromtxt
import pandas as pd
from osgeo import gdal
import os
#import xarray as xr
import clhs as cl
import csv
import numpy as np
from scipy.linalg import solve_triangular, get_lapack_funcs, get_blas_funcs
from maxvolpy.maxvol import maxvol


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


def f_no_cut(idx, i, copy=False):
    if copy:
        idx = np.copy(idx)
    idx[i] = 0
    return idx

def f_cut_eps(idx, i, X, eps=0.1, copy=False):
    if copy:
        idx = np.copy(idx)

    #print(np.abs(X - X[i]) < eps)
    print(idx.shape, X.shape)
    idx[np.abs(X - X[i]) < eps] = 0
    return idx

def rect_maxvol_cut(A, tol = 1., maxK = None, min_add_K = None, minK = None, start_maxvol_iters = 10, identity_submatrix = True, top_k_index = -1, cut_fun=None, penalty=None):
    """Python implementation of rectangular 2-volume maximization. For information see :py:func:`rect_maxvol` function"""
    # tol2 - square of parameter tol
    tol2 = tol**2
    # N - number of rows, r - number of columns of matrix A
    N, r = A.shape
    if N <= r:
        return np.arange(N, dtype = int), np.eye(N, dtype = A.dtype)
    if maxK is None or maxK > N:
        maxK = N
    if maxK < r:
        maxK = r
    if minK is None or minK < r:
        minK = r
    if minK > N:
        minK = N
    if min_add_K is not None:
        minK = max(minK, r + min_add_K) 
    if minK > maxK:
        minK = maxK
    if top_k_index == -1 or top_k_index > N:
        top_k_index = N
    if top_k_index < r:
        top_k_index = r

    if cut_fun is None:
        cut_fun = f_no_cut

    if penalty is None:
        #penalty_fun = np.ones(top_k_index, dtype=int)
        chosen = np.ones(top_k_index, dtype=int)
    else:
        chosen = np.copy(penalty)


    index = np.zeros(N, dtype = int)
    tmp_index, C = maxvol(A, tol = 1, max_iters = start_maxvol_iters, top_k_index = top_k_index)
    # -- 
    index[:r] = tmp_index
    #chosen[tmp_index] = 0 -- replaced
    for ti in tmp_index:
        cut_fun(chosen, ti)
    C = np.asfortranarray(C)
    # compute square 2-norms of each row in matrix C
    row_norm_sqr = np.array([chosen[i]*np.linalg.norm(C[i], 2)**2 for i in range(top_k_index)])
    # find maximum value in row_norm_sqr
    i = np.argmax(row_norm_sqr)
    K = r
    # set cgeru or zgeru for complex numbers and dger or sger for float numbers
    try:
        ger = get_blas_funcs('geru', [C])
    except:
        ger = get_blas_funcs('ger', [C])
    while (row_norm_sqr[i] > tol2 and K < maxK) or K < minK:
        # add i to index and recompute C and square norms of each row by SVM-formula
        index[K] = i
        #chosen[i] = 0 -- replaced by the next line
        #print(chosen)
        cut_fun(chosen, i)
        if (chosen == 0).all():
            print('Failed')
        c = C[i].copy()
        v = C.dot(c.conj())
        l = 1.0/(1+v[i])
        ger(-l,v,c,a=C,overwrite_a=1)
        C = np.hstack([C, l*v.reshape(-1,1)])
        row_norm_sqr -= (l*v[:top_k_index]*v[:top_k_index].conj()).real
        row_norm_sqr *= chosen
        # find maximum value in row_norm_sqr
        i = row_norm_sqr.argmax()
        K += 1
    if identity_submatrix:
        C[index[:K]] = np.eye(K, dtype = C.dtype)
    return index[:K].copy(), C


def make_dist(X):
    n = len(X)
    A = np.empty((n, n), dtype=X.dtype)
    for ix, x in enumerate(X):
        for iy, y in enumerate(X):
            A[ix, iy] = np.abs(x - y)

    return A

def f_penal(X, bnd, level=0.0):
    Xmin = np.min(X)
    Xmax = np.max(X)
    bnd_abs = (Xmax - Xmin)*bnd
    dist = np.minimum(np.abs(X - Xmin), np.abs(Xmax - X))
    def lin_func(x):
        if bnd == 0:
            return x*0.0 + 1.0  # crookedly, but it works. Ann, never do like this!
        else:
            return (1.0 - level)*np.minimum(x, bnd_abs)/bnd_abs + level

    return lin_func(dist)


def f_penal_2D(X, Y, bnd, level=0.0):
    return f_penal(X, bnd=bnd, level=level)*f_penal(Y, bnd=bnd, level=level)


def norm_data(X, bounds=(-1.0, 1.0), copy=True):
    X = np.array(X, copy=copy).T
    for i, x in enumerate(X):
        # print(len(x))
        min_v, max_v = np.min(x), np.max(x)
        b = (bounds[0]*max_v - bounds[1]*min_v)/(max_v-min_v)
        k = float(bounds[1] - bounds[0])/(max_v-min_v)
        X[i] *= k
        X[i] += b
        

    return X.T

def points_selection(X, max_n_pnts, min_n_pnts, cut_fun=None, penalty = None):
    
    """Function for selecting optimal parameters for dimentionality reduction method and for clustering.
    
    Parameters 
    ----------------
    X: array with shape (number_of_pixels*number_of_features)
            Initial data
           
    """
    #MaxVol
    
    res = rect_maxvol_cut(X, maxK=max_n_pnts, minK=min_n_pnts, cut_fun=cut_fun, penalty=penalty)[0]

    return res


def add_coords(X=None, size=(285, 217), order='C', idx_good_mask=None):
    """
    order can by 'C' or 'F'
    """
    w, h = size
    x_coord, y_coord = np.meshgrid(np.arange(h), np.arange(w))
    
    
    if idx_good_mask is None:
        idx_good_mask = np.arange(x_coord.size)
    
    if X is None:
        return np.hstack((
            x_coord.flatten(order=order)[idx_good_mask, None],
            y_coord.flatten(order=order)[idx_good_mask, None]))
    else:
        return np.hstack((np.array(X, copy=False),
                          x_coord.flatten(order=order)[idx_good_mask, None],
                          y_coord.flatten(order=order)[idx_good_mask, None]))
    
def gen_input(mode, data, shapes,mask):
    modes = ['usual', 'normed',
         'XY', 'XY_normed']
    fn_X_embedded = modes[mode]
    return [
        lambda x: np.array(x),
        lambda x: norm_data(x),
        lambda x: add_coords(
            x, size=shapes[0], idx_good_mask=mask),
        lambda x: norm_data(gen_input(2, x, shapes, mask)[0], copy=False),
    ][mode](data), fn_X_embedded

def my_score(a, b):
    a = np.array(a, copy=False)
    b = np.array(b, copy=False)
    n = len(a)
    assert len(b) == n, 'Arrays of differnet shapes :((('
    m = len(a[a==b])
    return float(m)/float(n)

def f_no_cut(idx, i, copy=False):
    if copy:
        idx = np.copy(idx)
    idx[i] = 0
    return idx

def f_cut_eps(idx, i, X, eps=0.1, copy=False):
    if copy:
        idx = np.copy(idx)
    xx = X[:, -2] 
    yy = X[:, -1]   
    #idx[i] = 0
    idx[(xx - xx[i])**2 + (yy-yy[i])**2 <= eps**2] = 0
    return idx


def calc_score(idx, X, y, to_ret_pred=True):
    gnb = GaussianNB()
    gnb_model = gnb.fit(X[idx], y[idx])
    
    if to_ret_pred:
        scores = extend_score(y, gnb_model.predict(X))
    else:
        scores = gnb_model.score(X, y)

    return scores


def good_points_brute_force(idx, num, X, y):
    sc = -1
    cmb_good = None
    for comb in combinations(idx, num):
        comb = np.array(comb)
        #print(comb)
        sc_curr = calc_score(comb, X=X, y=y, to_ret_pred=True)
        if sc_curr > sc:
            sc = sc_curr
            cmb_good = comb
            
    return cmb_good, sc

def idx_to_idx(idx_big, idx):
    hass = dict()
    for i, elem in enumerate(idx_big):
        hass[elem] = i
        
    return np.array([hass[i] for i in idx])



class MaxVolSampling():

    """
    Class to proccess data with MaxVol, cLHS and Random 
    
    Input: DEM, terrain features
    Return: Sampling points indices
    """

    def __init__(self):
        
        self.original_data = None
        self.maxvol_dist = None
        self.cLHS_dist = None
        self.random_dist = None
        self.maxvol_indices = None
        
        self.soil_feature = None
        self.num_of_points = 15
        self.path_to_file_with_indices = None
        self.wd = None
        self.soil_data = None
        self.X = None
        
        self.lons = None
        self.lats = None
        
        
        self.path_to_interpolation_file = None
        self.interpolation_array = None
    def data_preparation(self, wd, data_m, dem_dir):
        """
        Function to orginize tif files in flatten vectos, remove NaN and stack vectors into matrix

        """
        fl_names = list(filter(lambda fl: fl.endswith('.tif'), os.listdir(wd+'/ndvi_features/')))
        files = list(map(lambda x: gdal.Open(os.path.join(wd+'/ndvi_features/', x)), fl_names))
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

        dem_no_nodata = np.delete(dem_flat, idx_dem_nodata)

        #process with interp data
        if self.path_to_interpolation_file is not None:

            interpolation_raw_data = np.load(self.path_to_interpolation_file)[::-1]
            flatten_interpolation =  interpolation_raw_data.flatten()
            interpolation_no_nan = np.delete(flatten_interpolation, np.isnan(flatten_interpolation))
            self.interpolation_array = interpolation_no_nan

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
        self.X = X
    #     X = np.vstack((X, X[:1,:]))
        return X, dem_flat, dem_nodata, init_dem_shape, idx_dem, boundary_box


    def create_polygon(self, shape, vertices, value=1):
        """
        Creates np.array with dimensions defined by shape
        Fills polygon defined by vertices with ones, all other values zero"""
        base_array = np.zeros(shape, dtype=float)  # Initialize your array of zeros

        fill = np.ones(base_array.shape) * True  # Initialize boolean array defining shape fill

        # Create check array for each edge segment, combine into fill array
        for k in range(vertices.shape[0]):
            fill = np.all([fill, self.check(vertices[k-1], vertices[k], base_array)], axis=0)

        # Set all values inside polygon to one
        base_array[fill] = value

        return base_array,fill

    def find_nearest(self, array, value):
        array = np.asarray(array)
        idx = np.unravel_index(np.argmin((np.abs(array - value)), axis=None), array.shape)
        return array[idx], idx

    def check(self, p1, p2, base_array):
        """
        Uses the line defined by p1 and p2 to check array of 
        input indices against interpolated value

        Returns boolean array, with True inside and False outside of shape
        """
        idxs = np.indices(base_array.shape) # Create 3D array of indices

        p1 = p1.astype(float)
        p2 = p2.astype(float)

        # Calculate max column idx for each row idx based on interpolated line between two points

        if p1[0] == p2[0]:
            max_col_idx = (idxs[0] - p1[0]) * idxs.shape[1]
            sign = np.sign(p2[1] - p1[1])
        else:
            max_col_idx = (idxs[0] - p1[0]) / (p2[0] - p1[0]) * (p2[1] - p1[1]) + p1[1]
            sign = np.sign(p2[0] - p1[0])
        return idxs[1] * sign <= max_col_idx * sign  
    def original_soil_data(self, feature):
        soil_data = self.soil_data
        data = soil_data[feature]
        self.original_data = np.array(data)


    def dataframe_to_points(self):
        
        dem_raw = gdal.Open('dem.tif') #('/home/apetrovskaya/maxvol_soil_sampling/src/dem.tif')
        dem = dem_raw.ReadAsArray()
        self.init_dem_shape = dem.shape
        
        FEATURE = self.soil_feature
        
        soil_data = self.soil_data
        lons=soil_data['LON']
        self.lons = lons 
        lats=soil_data['LAT']
        self.lats = lats

        data = soil_data[FEATURE]
        self.original_data = np.array(data)
        #coordinate mesh
        xmin, ymin, xmax, ymax = [416949.0957, 5750852.2926,417891.8549,5751465.6945] #!!!HARDCODE
        st = dem
        xv = np.linspace(xmin,xmax, num=st.shape[1])
        yv = np.linspace(ymax,ymin, num=st.shape[0])
        coords = np.meshgrid(xv,yv)

        number_of_points=len(lons)
        points_idx=np.zeros((number_of_points,2))
        for i in range(number_of_points):
            a = self.find_nearest(coords[0],lons[i])[1][1]
            b = self.find_nearest(coords[1],lats[i])[1][0]
            points_idx[i,:]=[a,b]
            points_idx = points_idx.astype(int)
        return points_idx, data


    def distr_from_voronoi(self):

        points_idx,data = self.dataframe_to_points()

        #add points for right simplex
        points_idx_add = points_idx.copy()
        for i in range(-50,self.init_dem_shape[0]+50,50):
            points_idx_add = np.vstack((points_idx_add,[-50, i]))
            points_idx_add = np.vstack((points_idx_add,[self.init_dem_shape[1]+50,i]))

        for i in range(-50,self.init_dem_shape[1]+50,50):
            points_idx_add = np.vstack((points_idx_add,[i, -50]))
            points_idx_add = np.vstack((points_idx_add,[i,self.init_dem_shape[0]+50]))

        # generate Voronoi tessellation
        vor_add=Voronoi(points_idx_add)

        # cycle to fill regions in numpy array
        pol=np.zeros((self.init_dem_shape[1],self.init_dem_shape[0]))
        for r in range(len(vor_add.point_region)):
            region = vor_add.regions[vor_add.point_region[r]]

            if not -1 in region:
                value = data[r]
                polygon = [vor_add.vertices[i] for i in region]
                polygon = np.asarray(polygon)
                hull = ConvexHull(polygon)
                _, fill = self.create_polygon((self.init_dem_shape[1],self.init_dem_shape[0]),polygon[hull.vertices][::-1])
                pol[fill] = value
        pol[pol<min(data)]=min(data)
        polygons_in_array=pol.T
        self.voronoi_map = polygons_in_array.flatten()
        return self.voronoi_map
        
    
    
    def i_am_maxvol_function(self):
        
        self.num_of_points
        dist_pts = 0.1
        wd = self.wd
        data_m=3
        dem_dir = None
        max_n_pnts = self.num_of_points
        min_n_pnts = self.num_of_points

        X, dem_flat, dem_nodata, init_dem_shape, idx_dem, boundary_box = self.data_preparation(wd, data_m, dem_dir)

        #function for distance between points
        f_cut = lambda idx, i : f_cut_eps(idx, i, X=X, eps = dist_pts)
        #function for distance from border
        # f_penal = f_penal_2D(X = X[:, -2], Y = X[:, -1], bnd = 0.3, level = 0.3)
        f_penal = f_penal_2D(X = X[:, -2], Y = X[:, -1], bnd = 0.2, level = 0.3) #Change to 0.2 

        result = points_selection(X, max_n_pnts = max_n_pnts, min_n_pnts = min_n_pnts, cut_fun = f_cut, penalty = f_penal) 
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
        self.maxvol_indices = result
        return self.maxvol_indices


    def i_am_clhs(self, num_iter):

        n_pnts = self.num_of_points
        #cLHS
        sampled=cl.clhs(self.X[:,:-2], n_pnts, max_iterations=num_iter, progress=False)    
        self.cLHS_indices = sampled['sample_indices']
        return self.cLHS_indices 
    
    def i_am_random(self):
    
        random_dist = np.random.randint(low=0,high=self.X.shape[0],size=self.num_of_points)
        return random_dist        

if __name__ == "__main__":
    SAR = MaxVolSampling()
