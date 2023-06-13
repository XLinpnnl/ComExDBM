import csv
import sys
import pytz
from datetime import date, datetime, timedelta
from scipy import stats, signal, ndimage
import geopandas as gpd
import pyproj
from scipy.interpolate import griddata

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
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import dates
from matplotlib.font_manager import FontProperties
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.animation import FuncAnimation, PillowWriter
from colour import Color
import pylab as plt
import seaborn as sns
import imageio

import numpy.ma as ma
import xarray as xr
import rioxarray as rxr
from shapely.geometry import mapping, box
from shapely.geometry import Point, Polygon
import geopandas as gpd
import rasterio as rio
import warnings
#import metpy.calc as mpcalc
#from metpy.units import units
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

def MERRA2_asm(indir,outdir,year,region='central'):
    outfiles = [fn for fn in os.listdir(output_dir) if '.csv' in fn]
    ofile = 'MERRA_asm_'+str(year)+'.csv'
    #years= list(range(2001, 2021))
    #files= [fn for fn in os.listdir(indir) if ('.nc' in fn) and (str(year) in fn)]
    shapes = gpd.read_file(os.path.join('/qfs/projects/comexdbm/LDRD_extremes/ComExDBM/01_InputData','CONUS4326','contig_us.shp'), crs="EPSG:4326")
    states = shapes.NAME.unique()
    files= [fn for fn in os.listdir(indir) if ('.nc' in fn) and ('.'+str(year) in fn)]
    if ofile not in outfiles:
        fdata =[]
        for f in files:
            data = xr.open_dataset(os.path.join(indir,f))
            data.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
            data.rio.write_crs("EPSG:4326", inplace=True) #data.rio.set_crs("EPSG:4326") #
            #data['date']=data['time'].dt.strftime('%Y-%m-%d')#ds['Times'].str.decode("utf-8").str.replace('_',' ')#ds['Times'].str.slice(start=0, stop=10)
            #data['Td'] = mpcalc.dewpoint_from_relative_humidity(data['T'],data['RH']) #(data['T'] - 273.15)-((100-data['RH']*100)/5) #
            #data = data.groupby("date").mean();
            sset =[]
            for s in states:
                shape = shapes[shapes.NAME==s]
                sdata=  data.rio.clip(shape.geometry.apply(mapping), shape.crs, drop=True)
                isdf = sdata.to_dataframe().reset_index();isdf['date'] = pd.to_datetime(isdf['time']).dt.date
                isdf1 = isdf[isdf.lev==850].groupby(['date'],as_index=False).agg(
                        {'RH': ['mean','max'],
                         'T': ['mean','max'],
                        })	
                isdf1.columns=['date','RH_mean','RH_max','T_mean','T_max']
                isdf1['state'] = s.upper()
                ##
                isdf2= isdf[isdf.lev.isin([250,850])].groupby(['date','lev'],as_index=False).agg({'U': ['mean','max']})	
                isdf2.columns=['date','lev','U_mean','U_max']  
                isdf2 = isdf2.pivot(index='date', columns='lev', values=['U_mean','U_max']).reset_index()
                isdf2.columns=['date','U250_mean','U850_mean','U250_max','U850_max'] 
                isdf2['state'] = s.upper()
                idata =reduce(lambda left,right: pd.merge(left,right,on=['date','state'],how='outer'), [isdf1,isdf2])
                sset+=[idata]
            ##
            fdata+= [sset]
        ##
        fdata =pd.concat(fdata); 
        fdata = fdata.sort_values(by=['date'])
        #cdata =pd.concat(cdata)
        fdata.to_csv(os.path.join(output_dir,ofile),index=False)

##############################################################################
if __name__ == "__main__":
    #
    indir = '/qfs/projects/comexdbm/LDRD_extremes/'
    conus = gpd.read_file(os.path.join(indir,'CONUS4326','contig_us.shp'), crs="EPSG:4326")
    years= list(range(2001, 2021))
    months =[str(i).zfill(2) for i in list(range(1, 13))]
    #fdir=os.path.join(indir,'CPC_Global_Unified_Temperature')
    #tdata = xr.open_dataset(os.path.join(fdir,'tmax_tmin_tmean_conus_2001_2020.nc'))
    ndir=os.path.join(indir,'NARR_vars_daily')
    mdir=os.path.join(indir,'MERRA2_vars_daily')
    nfiles = [fn for fn in os.listdir(ndir)]
    mfiles = [fn for fn in os.listdir(mdir)]
    #regrid NARR data
    for y in years:
        #ydata = tdf[tdf.time.dt.year==y]
        outfile = 'NARR_vars_'+str(y)+'_regrid.csv'
        files= [fn for fn in os.listdir(os.path.join(ndir,'NARR_regrid'))]
        if outfile not in files:
            ndata = ndata = pd.read_parquet(os.path.join(ndir,'NARR_vars_'+str(y)+'.parquet'))
            ndata['date']=pd.to_datetime(ndata['date']);
            cols=ndata.columns.tolist()
            ydata = pd.read_csv(os.path.join(indir,'HW_Wildfire_Drought_data','CPC_Global_Unified_Temperature_with_HW_labels_Fires_Drought_'+str(y)+'.csv'))#.iloc[:, 1:]
            ydata['time']=pd.to_datetime(ydata['time']);
            lonlats = ydata.groupby(['lon', 'lat']).tmax.max().reset_index()
            ftdata = []
            for i in range(len(lonlats)):
                ilon = lonlats.lon.values[i]; ilat =lonlats.lat.values[i];
                idf = ydata[(ydata.lon==ilon)&(ydata.lat==ilat)]
                its = ndata[(ndata.lon>=(ilon-0.25))&(ndata.lon<=(ilon+0.25))&(ndata.lat>=(ilat-0.25))&(ndata.lat<=(ilat+0.25))]
                # its = its[its.date.dt.year==y]
                # its9 = its[its.conf==9]
                #its = its.groupby(['year']).gl.sum().reset_index()
                mcols1 = [c for c in ndata.columns if 'mean' in c]
                mcols2 = [c for c in ndata.columns if 'max' in c]
                mcols3 = [c for c in ndata.columns if 'min' in c]
                if (len(its)>0):
                    zs1 = its.groupby(['date'],as_index=True)[mcols1].mean().reset_index()
                    zs2 = its.groupby(['date'],as_index=True)[mcols2].max().reset_index()
                    zs3 = its.groupby(['date'],as_index=True)[mcols3].min().reset_index()
                    zs = reduce(lambda x, y: pd.merge(x, y, on = 'date'), [zs1,zs2,zs3])
                    zs['lon'] =ilon; zs['lat'] =ilat; 
                    zs=zs[cols]
                    ftdata +=[zs]
            #
            ftdata=pd.concat(ftdata)  
            ftdata.to_csv(os.path.join(ndir,'NARR_regrid','NARR_vars_'+str(y)+'.csv'),index=False)

    #regrid MERRA data
    for y in years:
        #ydata = tdf[tdf.time.dt.year==y]
        outfile = 'MERRA_asm_'+str(y)+'_regrid.csv'
        files= [fn for fn in os.listdir(os.path.join(mdir,'MERRA_regrid'))]
        if outfile not in files:
            ndata =  pd.read_parquet(os.path.join(mdir,'MERRA_asm_'+str(y)+'.parquet'))
            ndata['date']=pd.to_datetime(ndata['date']);
            cols=ndata.columns.tolist()
            ydata = pd.read_csv(os.path.join(indir,'HW_Wildfire_Drought_data','CPC_Global_Unified_Temperature_with_HW_labels_Fires_Drought_'+str(y)+'.csv'))#.iloc[:, 1:]
            ydata['time']=pd.to_datetime(ydata['time']);
            lonlats = ydata.groupby(['lon', 'lat']).tmax.max().reset_index()
            ftdata = []
            for i in range(len(lonlats)):
                ilon = lonlats.lon.values[i]; ilat =lonlats.lat.values[i];
                idf = ydata[(ydata.lon==ilon)&(ydata.lat==ilat)]
                its = ndata[(ndata.lon>=(ilon-0.25))&(ndata.lon<=(ilon+0.25))&(ndata.lat>=(ilat-0.25))&(ndata.lat<=(ilat+0.25))]
                # its = its[its.date.dt.year==y]
                # its9 = its[its.conf==9]
                #its = its.groupby(['year']).gl.sum().reset_index()
                mcols1 = [c for c in ndata.columns if 'mean' in c]
                mcols2 = [c for c in ndata.columns if 'max' in c]
                mcols3 = [c for c in ndata.columns if 'min' in c]
                if (len(its)>0):
                    zs1 = its.groupby(['date'],as_index=True)[mcols1].mean().reset_index()
                    zs2 = its.groupby(['date'],as_index=True)[mcols2].max().reset_index()
                    zs3 = its.groupby(['date'],as_index=True)[mcols3].min().reset_index()
                    zs = reduce(lambda x, y: pd.merge(x, y, on = 'date'), [zs1,zs2,zs3])
                    zs['lon'] =ilon; zs['lat'] =ilat; 
                    zs=zs[cols]; zs.columns = [c.replace('_cow','') for c in zs.columns]
                    ftdata +=[zs]
            #
            ftdata=pd.concat(ftdata)  
            ftdata.to_csv(os.path.join(mdir,'MERRA_regrid','MERRA_asm_'+str(y)+'.csv'),index=False)
    
    # ndir=os.path.join(indir,'NARR_vars_daily','NARR_regrid')
    # mdir=os.path.join(indir,'MERRA2_vars_daily','MERRA_regrid')
    # nfiles = [fn for fn in os.listdir(ndir)]
    # mfiles = [fn for fn in os.listdir(mdir)]
    # for f in nfiles:
        # ndata =  pd.read_csv(os.path.join(ndir,f))
        # ndata.to_parquet(os.path.join(ndir,f.replace('.csv','.parquet')),index=False)

    # #
    # for f in mfiles:
        # mdata =  pd.read_csv(os.path.join(mdir,f))
        # mdata.to_parquet(os.path.join(mdir,f.replace('.csv','.parquet')),index=False)
    # for f in range(len(nfiles)):
        # ndata = pd.read_parquet(os.path.join(ndir,nfiles[f]))
        # ndata['date'] = pd.to_datetime(ndata['date'])
        # mcols= [c for c in gdf.columns if 'mean' in c]
        # df_rows = pd.DataFrame(ndata).set_index(["lon", "lat","date"])
        # df_rows = df_rows[mcols]
        # data = xr.Dataset.from_dataframe(df_rows)
        # data.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
        # data.rio.write_crs("EPSG:4326", inplace=True)
        # #data=  data.rio.clip(conus.geometry.apply(mapping), conus.crs, drop=True)
        # #data = xr.open_dataset(os.path.join(indir,index+'.nc'))
        # #data = data.sel(day=slice('2001-01-01', '2020-12-31'))
        # data_out = xr.Dataset(
            # {
                # "lat": (["lat"], np.arange(25.25, 49.25, 0.5)),
                # "lon": (["lon"], np.arange(-124.75, -66.75, 0.5)),
            # }
        # )
        # regridder = xe.Regridder(data, data_out, "bilinear")
        # data_out = regridder(data)
        # #data_out.attrs=data.attrs
        # data_out.to_netcdf(os.path.join(ndir,'NARR_regrid',nfiles[f].replace('.parquet','_mean.nc'))
    # ###
    # for f in range(len(mfiles)):
        # ndata = pd.read_parquet(os.path.join(ndir,mfiles[f]))
        # ndata['date'] = pd.to_datetime(ndata['date'])
        # mcols= [c for c in gdf.columns if 'mean' in c]
        # df_rows = pd.DataFrame(ndata).set_index(["lon", "lat","date"])
        # df_rows = df_rows[mcols]
        # data = xr.Dataset.from_dataframe(df_rows)
        # data.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
        # data.rio.write_crs("EPSG:4326", inplace=True)
        # #data=  data.rio.clip(conus.geometry.apply(mapping), conus.crs, drop=True)
        # #data = xr.open_dataset(os.path.join(indir,index+'.nc'))
        # #data = data.sel(day=slice('2001-01-01', '2020-12-31'))
        # data_out = xr.Dataset(
            # {
                # "lat": (["lat"], np.arange(25.25, 49.25, 0.5)),
                # "lon": (["lon"], np.arange(-124.75, -66.75, 0.5)),
            # }
        # )
        # regridder = xe.Regridder(data, data_out, "bilinear")
        # data_out = regridder(data)
        # #data_out.attrs=data.attrs
        # data_out.to_netcdf(os.path.join(mdir,'MERRA_regrid',nfiles[f].replace('.parquet','_mean.nc'))
