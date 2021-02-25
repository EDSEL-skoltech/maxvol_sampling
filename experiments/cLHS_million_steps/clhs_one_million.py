from scipy.spatial import voronoi_plot_2d, Voronoi
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull, convex_hull_plot_2d
from numpy import genfromtxt
import pandas as pd
import numpy as np
import gdal
import argparse


import os
import xarray as xr
#cLHS
import clhs as cl
import csv


from scipy.spatial import distance
from scipy.stats import entropy
from scipy.special import kl_div



import argparse
import numpy as np
import osgeo.gdal as gdal
import os
#by Anna
from data_preparation import data_preparation
from maxvol_cut import rect_maxvol_cut, f_no_cut, f_penal_2D
from tools import norm_data, add_coords, gen_input, extend_score, points_selection, f_no_cut, f_cut_eps, calc_score, good_points_brute_force, idx_to_idx
import csv


class Divergence():

    """
    class to proccess data with MaxVol, cLHS and Random 
    
    I hope it will return dist
    """

    def __init__(self):
        
        self.original_data = None
        self.maxvol_dist = None
        self.cLHS_dist = None
        self.random_dist = None
        self.maxvol_indices = None
        self.cLHS_indices = None
        
        self.soil_feature = None
        self.num_of_points = 15
        self.path_to_file_with_indices = None
        self.wd = None
        self.soil_data = None
        self.X = None
        
    def data_preparation(self, wd, data_m, dem_dir):
        """
        Function to orginize tif files in flatten vectos, remove NaN and stack vectors into matrix

        """
        fl_names = list(filter(lambda fl: fl.endswith('.tif'), os.listdir(wd+'/features/')))
        files = list(map(lambda x: gdal.Open(os.path.join(wd+'/features/', x)), fl_names))
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
        self.X = X
    #     X = np.vstack((X, X[:1,:]))
        return X, dem_flat, dem_nodata, init_dem_shape, idx_dem


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

    def dataframe_to_points(self):
        
        dem_raw = gdal.Open('dem.tif') #('/home/apetrovskaya/maxvol_soil_sampling/src/dem.tif')
        dem = dem_raw.ReadAsArray()
        self.init_dem_shape = dem.shape
        
        FEATURE = self.soil_feature
        
        soil_data = self.soil_data
        
        
        lons=soil_data['LON']
        lats=soil_data['LAT']
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


    def distr_from_voronoi(self, sampling_result):

#         points_idx,data = self.dataframe_to_points(dataframe, feature)
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

        #delete 0s
        pol[pol<min(data)]=min(data)
        polygons_in_array=pol.T

        distribution = polygons_in_array.flatten()[sampling_result]

        return distribution
    

    def i_am_maxvol_function(self):
        
        self.num_of_points
        dist_pts = 0.3

        wd = '../src/'
        data_m=3
        dem_dir = None

        max_n_pnts = self.num_of_points

        min_n_pnts = self.num_of_points

        X, dem_flat, dem_nodata, init_dem_shape, idx_dem = data_preparation(wd, data_m, dem_dir)

        print('Selecting points...')
        #function for distance between points
        f_cut = lambda idx, i : f_cut_eps(idx, i, X=X, eps = dist_pts)
        #function for distance from border
        f_penal = f_penal_2D(X = X[:, -2], Y = X[:, -1], bnd = 0.3, level = 0.3)

        result = points_selection(X, max_n_pnts = max_n_pnts, min_n_pnts = min_n_pnts, cut_fun = f_cut, penalty = f_penal) 

        #coordinates
        xmin, ymin, xmax, ymax = [37.7928399,51.90236556, 37.8064010,51.90774268]

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

        # np.savetxt(args.wd+'/maxvol_result_'+str(len(result))+'.csv', coord_idx, delimiter=';')
        # np.save(args.wd+'/maxvol_idx_'+str(len(result))+'.npy', result)

        # ## add by Misha
        # with open('temp_result_for_original.csv', 'a') as f:
        #     writer = csv.writer(f)
        #     writer.writerow(result)
        
        self.maxvol_indices = result
        
        self.maxvol_dist = self.distr_from_voronoi(result)
        return self.maxvol_dist


    def i_am_clhs(self, num_iter):
        wd = '../src/'
        dem_dir = None
        n_pnts = self.num_of_points

        #cLHS
        sampled=cl.clhs(self.X[:,:4], n_pnts, max_iterations=num_iter)
        self.cLHS_indices = sampled['sample_indices']
        
        self.cLHS_dist = self.distr_from_voronoi(sampled['sample_indices'])

        return self.cLHS_dist 
    
    def i_am_random(self):
    
        
        random_dist = np.random.randint(low=0,high=self.X.shape[0],size=self.num_of_points)
        self.random_dist = self.distr_from_voronoi(random_dist)
    
        return self.random_dist
if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--num_points', type=int, help='num of points', required=True)
    args = parser.parse_args()

    SAR = Divergence()
    SAR.soil_feature = 'Moisture_perc_30'
    SAR.num_of_points = args.num_points

    SAR.soil_data = pd.read_csv('kshen_new_data.csv', sep = ";")

    SAR.path_to_file_with_indices = None
    SAR.wd = '../src'

    _  =SAR.data_preparation(SAR.wd, data_m=3, dem_dir = None)

    cLHS = SAR.i_am_clhs(num_iter=1000000)

    with open('./one_million_cLHS'+str(args.num_points)+'_5m_DEM.csv', 'w') as f:
        writer = csv.writer(f)
        writer.writerow(SAR.cLHS_indices)    



    print("____DONE_____")