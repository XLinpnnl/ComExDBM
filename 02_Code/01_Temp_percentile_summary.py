import csv
import sys
import pytz
from datetime import date, datetime, timedelta
from scipy import stats, signal, ndimage

import os
import re
import netCDF4
import numpy as np
import pandas as pd
from functools import reduce

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.text as text
import matplotlib.ticker as mticker
import matplotlib.colors as mcolors
from matplotlib import dates
from matplotlib.font_manager import FontProperties
import pylab as plt
import seaborn as sns

import numpy.ma as ma
import xarray as xr
import rioxarray as rxr
from shapely.geometry import mapping, box
from shapely.geometry import Point, Polygon
import geopandas as gpd
import rasterio as rio
import warnings
import metpy.calc as mpcalc
from metpy.units import units
import itertools
import warnings
warnings.filterwarnings("ignore")

#================================================================================
#shp_path = os.path.join('./','great_plains/great_plains.shp')

def fminus(x):
	return min(x)-min(x)

# lons = data.lon.unique();lats = data.lat.unique()
# lon = lons[0];lat = lats[0]
# data[(data.lon==lon)&(data.lon==lon)]
def sub_index(datatime,end):
	start =end- np.timedelta64(30,'D')
	return np.where((datatime>=start)&(datatime<=end))[0]
#ind_list = list(map(lambda i:sub_index(datatime, i), t0))

def dist(lat1, long1, lat2, long2):
	"""
	Calculate the great circle distance between two points 
	on the earth (specified in decimal degrees)
	"""
	# convert decimal degrees to radians 
	lat1, long1, lat2, long2 = map(radians, [lat1, long1, lat2, long2])
	# haversine formula 
	dlon = long2 - long1 
	dlat = lat2 - lat1 
	a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
	c = 2 * asin(sqrt(a)) 
	# Radius of earth in kilometers is 6371
	km = 6371* c
	return km

def find_nearest(lat, long):
	distances = edata.apply(lambda row: dist(lat, long, row['lat'], row['lon']), axis=1)
	return edata.loc[distances.idxmin(), ['slope', 'northness', 'elevation']]

def corners_to_poly(hdf):
	"""Convert daily haildata to polygon
	"""
	min_y= np.min([hdf.BEGIN_LAT.min(),hdf.END_LAT.min()])
	min_y = np.min([hdf.BEGIN_LAT.min(),hdf.END_LAT.min()])
	min_x= np.min([hdf.BEGIN_LON.min(),hdf.END_LON.min()])
	min_x = np.min([hdf.BEGIN_LON.min(),hdf.END_LON.min()])
	# create shapely polygon from corner coordinates
	geom = Polygon(
		((min_x, min_y),
		 (min_x, min_y),
		 (min_x, min_y),
		 (min_x, min_y),
		 (min_x, min_y))
		 )
	return geom

# def point_in_state(shape,state,pdata):
    # indx =[n for n in range(len(shape.NAME)) if shape.NAME[n] in state]
    # shape = shape.iloc[indx]
    # #point_all = Point(pdata.lon.values, pdata.lat.values)
    # cols = pdata.columns.tolist()
    # points = gpd.GeoDataFrame(pdata, geometry=gpd.points_from_xy(pdata.lon, pdata.lat))
    # points.crs = {'init': 'epsg:4326'}
    # within_points = gpd.sjoin(points, shape, op = 'within')
    # wdata = within_points[cols]
    # #wdata = pd.DataFrame(wdata.drop('geometry',axis=1))
    # #wdata.to_parquet(os.path.join(path,'precip_monthly_indices_conus','pmin3h_'+region+'_parquet_updated0301.gzip'), compression='gzip', engine='pyarrow')
    # return wdata 
def point_in_state(shapes,state,pdata):
    #indx =[n for n in range(len(shapes.NAME)) if shapes.NAME[n] in state]
    #shape = shape.iloc[indx]
    shape =shapes[shapes.NAME==state]
    #point_all = Point(pdata.lon.values, pdata.lat.values)
    cols = pdata.columns.tolist()
    points = gpd.GeoDataFrame(pdata, geometry=gpd.points_from_xy(pdata.lon, pdata.lat))
    points.crs = {'init': 'epsg:4326'}
    within_points = gpd.sjoin(points, shape, op = 'within')
    wdata = within_points[cols]
    wdata['state'] =state
    #wdata = pd.DataFrame(wdata.drop('geometry',axis=1))
    #wdata.to_parquet(os.path.join(path,'precip_monthly_indices_conus','pmin3h_'+region+'_parquet_updated0301.gzip'), compression='gzip', engine='pyarrow')
    return wdata 

def GUTemp_state_daily_tmax_tmin(tdir,year,state,shapes):
    # read in CPC Global Unified Temperature(e.g.,tmin,tmax) by state
    # calculate daily men (tmean = (tmax+tmin)/2)
    ftmins = [fn for fn in os.listdir(tdir) if 'tmin' in fn]
    ftmins.sort()
    ftmaxs = [fn for fn in os.listdir(tdir) if 'tmax' in fn]
    ftmaxs.sort()
    shape = shapes[shapes.NAME==state]
    bounds = shape.geometry.bounds
    tmax = xr.open_dataset(os.path.join(tdir,[f for f in ftmaxs if str(year) in f][0]))
    # tmin.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
    # tmin.rio.write_crs("EPSG:4326", inplace=True) #data.rio.set_crs("EPSG:4326") 
    tmin = xr.open_dataset(os.path.join(tdir,[f for f in ftmins if str(year) in f][0]))
    # tmin.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
    # tmin.rio.write_crs("EPSG:4326", inplace=True) #data.rio.set_crs("EPSG:4326") 
    data = xr.merge([tmax,tmin])
    data['tmean'] = (data['tmax']+data['tmin'])/2
    data = data.assign_coords(lon=(((data.lon + 180) % 360) - 180))
    #data.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
    sdata = data.where((data.lon>=bounds.minx.values)&(data.lon<=bounds.maxx.values)&(data.lat>=bounds.miny.values)&(data.lat<=bounds.maxy.values), drop=True) #-125,25,-66,50
    sdata.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
    sdata.rio.write_crs("EPSG:4326", inplace=True)
    sdata0=  sdata.rio.clip(shape.geometry.apply(mapping), shape.crs, drop=True)
    #tmin['Date'] = str(y)+'-'+tmin.Month.apply(str).str.zfill(2)+'-'+tmin.Day.apply(str).str.zfill(2)
    isdf = sdata.to_dataframe().reset_index();    
    isdf = isdf.dropna(); isdf=isdf.drop('spatial_ref',axis=1)
    isdf = isdf.rename(columns={'time': 'date'})
    isdf['state'] =state
    return isdf


##############################################################################
if __name__ == "__main__":
    ##PNW(WA,ID,OR)
    # file system
    indir = '/home/linx882/LDRD_extremes/'
    conus = gpd.read_file(os.path.join(indir,'CONUS4326','contig_us.shp'), crs="EPSG:4326")
    years= list(range(2001, 2021))
    PNW_shape = conus[conus.NAME.isin(['Washington','Idaho','Oregon'])]

    #########################################################
    ## CPC Global Unified Temperature
    tdir = os.path.join(indir,'CPC_Global_Unified_Temperature') #os.getcwd()
    ftmins = [fn for fn in os.listdir(tdir) if 'tmin' in fn]
    ftmins.sort()
    ftmaxs = [fn for fn in os.listdir(tdir) if 'tmax' in fn]
    ftmaxs.sort()
    PNW_shape = conus[conus.NAME.isin(['Washington','Idaho','Oregon','California'])]
    states = conus.NAME.unique().tolist()
    PNW= ['Washington','Idaho','Oregon','California']

    ## get temparature percentiles for each grid
    cdata =[]
    for y in years:
        dmin = xr.open_dataset(os.path.join(tdir,[f for f in ftmins if str(y) in f][0]))
        dmax = xr.open_dataset(os.path.join(tdir,[f for f in ftmaxs if str(y) in f][0]))
        data = xr.merge([dmax,dmin])
        data['tmean'] = (data['tmax']+data['tmin'])/2
        # tmin.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
        # tmin.rio.write_crs("EPSG:4326", inplace=True) #data.rio.set_crs("EPSG:4326") 
        #tmin = xr.open_dataset(os.path.join(tdir,[f for f in ftmins if str(year) in f][0]))
        # tmin.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
        # tmin.rio.write_crs("EPSG:4326", inplace=True) #data.rio.set_crs("EPSG:4326") 
        data = data.assign_coords(lon=(((data.lon + 180) % 360) - 180))
        #data.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
        sdata = data.where((data.lon>=-125)&(data.lon<=-66)&(data.lat>=25)&(data.lat<=50), drop=True) #-125,25,-66,50
        sdata.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
        sdata.rio.write_crs("EPSG:4326", inplace=True)
        sdata=  sdata.rio.clip(conus.geometry.apply(mapping), conus.crs, drop=True)
        cdata +=[sdata]

    kcdata= xr.concat(cdata,dim ='time')
    kcdata.to_netcdf(os.path.join(tdir,'tmax_tmin_tmean_conus_2001_2020.nc'))
    #
    kcdata = xr.open_dataset(os.path.join(tdir,'tmax_tmin_tmean_conus_2001_2020.nc'))
    qt_dims = ("time")
    qt_values = (0.75,0.81, 0.9,0.95,0.975, 0.98, 0.99)
    isdf = kcdata.quantile(qt_values, dim=qt_dims)
    #tmin['Date'] = str(y)+'-'+tmin.Month.apply(str).str.zfill(2)+'-'+tmin.Day.apply(str).str.zfill(2)
    isdf = isdf.to_dataframe().reset_index(); isdf = isdf.dropna(); 
    #isdf=isdf.drop('spatial_ref',axis=1)
    #isdf['quantile'] = isdf['quantile']*100#.astype(int)
    data = isdf.pivot(index=['lon','lat'], columns='quantile', values=['tmin','tmax','tmean'])
    data.columns=[s1 + str(s2*100).replace('.0', '') for (s1,s2) in data.columns.tolist()]
    data= data.reset_index()
    data.to_csv(os.path.join(tdir,'Temp_thresholds.csv'))


# #########################################################
# #MODIS fire data
# cwd =  os.path.join(indir,'MOD14A1') #os.getcwd()
# cp = pd.read_hdf(os.path.join(cwd, 'cp.h5'))
# # load fire events and component information
# #v = pd.read_hdf(os.path.join(cwd, 'v.h5'))
# v = pd.read_hdf(
    # os.path.join(cwd, 'v.h5'),
# )
# cols = ['t', 'dtime', 'lat', 'lon', 'area', 'maxFRP', 'neigh_int', 'gl', 'cp','conf','neigh']
# points = gpd.GeoDataFrame(v[cols], geometry=gpd.points_from_xy(v.lon, v.lat))
# points.crs = {'init': 'epsg:4326'}
# vps = gpd.sjoin(points, conus, op = 'within')
# vps = vps [cols+['NAME']]
# vps = vps.rename(columns={'NAME':'state'})

# vp =vps[vps.conf==9]
# # v=v.drop('state',axis=1)
# # v = v.merge(cp[['cp','state']],on='cp',how='left')
# # v= v.dropna(subset=['state'],axis=0)

# PNW= ['Washington','Idaho','Oregon']
# sv3 = vp[vp.state.isin(PNW)]
# sv3['dtime'] = pd.to_datetime(sv3['dtime'], format='%Y-%m-%d') #format='%Y-%m-%d' '%m/%d/%Y'
# sv3['year'] =sv3['dtime'].dt.year; sv3['year']=sv3['year'].astype('int')
# sgls = sv3.groupby(['lat', 'lon']).gl.count().reset_index()

# yts = sv3.groupby(['lat', 'lon','year']).gl.count().reset_index()
# yts[(yts.gl>5)&(yts.gl<10)]

# #data=sgls
# for year in  years:
    # data = yts[yts.year==year]
    # fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(6, 6))
    # #ax = conus_new.plot(facecolor="none",edgecolor=conus['region'],linewidth=2)
    # ax.set_title('Number of fire detection in each grid',fontsize=8)
    # divider = make_axes_locatable(ax)
    # cax = divider.append_axes("right", size="3%", pad=0.01)
    # x, y = data['lon'].values, data['lat'].values
    # # x, y = tdata_ngp['lon'].values, tdata_ngp['lat'].values
    # im=ax.scatter(x, y, c=data['gl'], marker="s",s=1,cmap=plt.get_cmap("coolwarm"),alpha=0.6) #,marker="s",alpha=0.8, cmap=color_ramp
    # #im=ax.scatter(x, y, marker="*",s=0.01) #,marker="s"
    # #im.set_clim(data['gl'].quantile(.1),data['gl'].quantile(.9))
    # im.set_clim(1,3)
    # #conus.plot(facecolor="none",edgecolor='r',linewidth=1,alpha=0.2,ax=ax)
    # PNW_shape.plot(facecolor="none",edgecolor='k',linewidth=1,ax=ax)
    # fig.colorbar(im, cax=cax, orientation='vertical',shrink=0.9)
    # ax.set_xlabel('Lon')
    # ax.set_ylabel('Lat')
    # plt.savefig(os.path.join('/home/linx882/LDRD_extremes/','fire_plots','WildFire_summary_grids_'+str(year)+'.png'), bbox_inches='tight',dpi=300)






