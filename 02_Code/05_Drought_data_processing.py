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
import metpy.calc as mpcalc
from metpy.units import units
import itertools
import warnings
import xarray as xr
import xesmf as xe
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

def regrid_drought_index(indir, index):
    data = xr.open_dataset(os.path.join(indir,index+'.nc'))
    data = data.sel(day=slice('2001-01-01', '2020-12-31'))
    data_out = xr.Dataset(
        {
            "lat": (["lat"], np.arange(25.25, 49.25, 0.5)),
            "lon": (["lon"], np.arange(-124.75, -67.25, 0.5)),
        }
    )
    regridder = xe.Regridder(data, data_out, "bilinear")
    data_out = regridder(data)
    data_out.attrs=data.attrs
    data_out.to_netcdf(os.path.join(ddir,index+'_regrid.nc'))
    return data_out
  
##############################################################################
if __name__ == "__main__":
    ##PNW(WA,ID,OR)
    # file system
    indir = '/home/linx882/LDRD_extremes/'
    conus = gpd.read_file(os.path.join(indir,'CONUS4326','contig_us.shp'), crs="EPSG:4326")
    years= list(range(2001, 2021))
    months =[str(i).zfill(2) for i in list(range(1, 13))]
    PNW_shape = conus[conus.NAME.isin(['Washington','Idaho','Oregon','California'])]
    #########################################################
    ## CPC Global Unified Temperature
    tdir = os.path.join(indir,'CPC_Global_Unified_Temperature') #os.getcwd()
    states = conus.NAME.unique().tolist()
    PNW= ['Washington','Idaho','Oregon','California']
    #  #################################################################################################
    # regrid drought data to 0.5°×0.5° spatial resolution
    ddir = os.path.join(indir,'Drought_data/gridMET')
    dindex = ['pdsi','spi14d','spi30d','spei14d','spei30d']
    for d in dindex :
        idata = regrid_drought_index(indir=ddir, index=d)
    
    
    # assign category to drought indices 
    years= list(range(2001, 2021))
    ydf = []
    for y in years:
        #ydata = tdf[tdf.time.dt.year==y]
        #fdata['discovery_date']=pd.to_datetime(fdata['discovery_date']);
        yfile = 'CPC_Global_Unified_Temperature_with_HW_labels_Fires_'+str(y)+'.csv'
        ofile = 'CPC_Global_Unified_Temperature_with_HW_labels_Fires_Drought_'+str(y)+'.csv'
        ydata = pd.read_csv(os.path.join(indir,'HW_Wildfire_Drought_data',yfile))#.iloc[:, 1:]
        ydata['time']=pd.to_datetime(ydata['time']);
        dcol = [c for c in ydata.columns if 'Unnamed' in c]
        if len(dcol)>0:
            ydata = ydata.drop(dcol,axis=1)
        #
        dindex = ['pdsi','spi14d','spi30d','spi90d','spei14d','spei30d','spei90d']
        for idx in dindex :
            #data = xr.open_dataset(os.path.join(indir,'Drought_data/gridMET',idx+'.nc'))
            data = xr.open_dataset(os.path.join(indir,'Drought_data/gridMET',idx+'_regrid.nc'))
            idata = data.sel(day=slice(str(y)+'-01-01', str(y)+'-12-31'))
            idf = idata.to_dataframe().reset_index(); idf = idf.dropna(); 
            idf.category = idf.category.astype(int)
            svar = idf.columns[4]
            if svar=='pdsi':
                idf.category[idf[svar]<= -5] = 0
                idf.category[(idf[svar]>-5)&(idf[svar]<=-4)] = 1
                idf.category[(idf[svar]>-4)&(idf[svar]<=-3)] = 2
                idf.category[(idf[svar]>-3)&(idf[svar]<=-2)] = 3
                idf.category[(idf[svar]>-2)&(idf[svar]<=-1)] = 4
                idf.category[(idf[svar]>-1)&(idf[svar]<1)] = 5
                idf.category[(idf[svar]>=1)&(idf[svar]<2)] = 6
                idf.category[(idf[svar]>=2)&(idf[svar]<3)] = 7
                idf.category[(idf[svar]>=3)&(idf[svar]<4)] = 8
                idf.category[(idf[svar]>=4)&(idf[svar]<5)] = 9
                idf.category[(idf[svar]>=5)] = 10
            else:
                idf.category[idf[svar]<= -2] = 0
                idf.category[(idf[svar]>-2)&(idf[svar]<=-1.6)] = 1
                idf.category[(idf[svar]>-1.6)&(idf[svar]<=-1.3)] = 2
                idf.category[(idf[svar]>-1.3)&(idf[svar]<=-0.8)] = 3
                idf.category[(idf[svar]>-0.8)&(idf[svar]<=-0.5)] = 4
                idf.category[(idf[svar]>-0.5)&(idf[svar]<0.5)] = 5
                idf.category[(idf[svar]>=0.5)&(idf[svar]<0.8)] = 6
                idf.category[(idf[svar]>=0.8)&(idf[svar]<1.3)] = 7
                idf.category[(idf[svar]>=1.3)&(idf[svar]<1.6)] = 8
                idf.category[(idf[svar]>=1.6)&(idf[svar]<2)] = 9
                idf.category[(idf[svar]>=2)] = 10
            ##
            idf.columns = ['time', 'lat', 'lon', 'crs', idx, idx+'_category']
            idf= idf.drop(['crs'],axis=1)
            start = datetime.strptime(str(y)+'-01-01', "%Y-%m-%d")
            end = datetime.strptime(str(y+1)+'-01-01', "%Y-%m-%d")
            date = pd.date_range(start,end-timedelta(days=1),freq='d')
            lonlat = idf.groupby(['lat', 'lon']).time.min().reset_index()
            #nll= lonlat.loc[lonlat.index.repeat(len(date))].reset_index(drop=True)
            nnl = pd.concat([lonlat]*len(date))
            ndate = date.repeat(len(lonlat)); nnl.time = ndate
            ndf =nnl.merge(idf,on =['time', 'lat', 'lon'],how='left')
            ndf[[idx, idx+'_category']] =ndf.groupby(['lat', 'lon'])[[idx, idx+'_category']].fillna(method='bfill')
            ydata = ydata.merge(ndf,on =['lat', 'lon','time'],how='left')
        #
        ydata.to_csv(os.path.join(indir,'HW_Wildfire_Drought_data',ofile),index=False)
    # #pdata =pd.read_csv(os.path.join(tdir,'Tmax_Tmin_thresholds.csv'))
    # pdata= pd.read_csv(os.path.join(tdir,'Temp_thresholds.csv'))
    # tdata = xr.open_dataset(os.path.join(tdir,'tmax_tmin_tmean_conus_2001_2020.nc'))
    # tdf = tdata.to_dataframe().reset_index(); #tdf = tdf.dropna(); 
    # tdf =tdf.drop('spatial_ref',axis=1)
    # tdf =tdf.merge(pdata[['lon', 'lat','state']],on=['lon', 'lat'],how='left')
    # tdf=tdf[tdf.state.notna()]

    #  #################################################################################################
    # years= list(range(2011, 2021))
    # fod = gpd.read_file(os.path.join(indir,'FPA_FOD_20210617','FPA_FOD_20210617.gpkg'))
    # fod.columns = [x.lower() for x in fod.columns]
    # fod=fod.rename({'discovery_date': 'disc_date','latitude':'lat','longitude':'lon'}, axis=1)
    # fod['disc_date']=fod.disc_date.apply(lambda x: x.date());#fod['cont_date']=fod.cont_date.apply(lambda x: x.date()) 
    # #fod['cont_date']=pd.to_datetime(fod['cont_date'],format='%Y-%m-%dT%H:%M:%S')
    # fod['disc_date']=pd.to_datetime(fod['disc_date']);
    # fod=fod[fod.disc_date.dt.year>=2001]
    # #fod.cont_date[fod.fod_id==1336443] = '2001-06-19T00:00:00+00:00'
    # fod['end_date']=pd.to_datetime(fod['cont_date'].str[:10], errors = 'coerce');
    # year = fod['disc_date'].dt.year.astype(str)
    # month_day = fod['end_date'].dt.strftime('%m-%d')
    # fod['end_date'][(fod.end_date-fod.disc_date)/np.timedelta64(1, 'D')>30] = pd.to_datetime(year +'-'+ month_day, format='%Y-%m-%d')[(fod.end_date-fod.disc_date)/np.timedelta64(1, 'D')>30]
    # fod['end_date'][fod['end_date'].isna()] = fod['disc_date'][fod['end_date'].isna()]
    # fod.fire_size = fod.fire_size/247.11
    # fod = fod[['fod_id', 'fpa_id','fire_year', 'disc_date','end_date','cont_date','fire_size', 'fire_size_class', 'lat', 'lon']]
    # for y in years:
        # #ydata = tdf[tdf.time.dt.year==y]
        # fdata = fod[fod.fire_year==y]
        # #fdata['discovery_date']=pd.to_datetime(fdata['discovery_date']);
        # outfile = 'CPC_Global_Unified_Temperature_with_HW_labels_Fires_'+str(y)+'.csv'
        # ydata = pd.read_csv(os.path.join(indir,'HW_Wildfire_Drought_data','CPC_Global_Unified_Temperature_with_HW_labels_Fires_'+str(y)+'.csv'))#.iloc[:, 1:]
        # ydata= ydata.drop(['fire_size_FOD'],axis=1)
        # ydata['time']=pd.to_datetime(ydata['time']);
        # lonlats = ydata.groupby(['lon', 'lat']).tmax.max().reset_index()
        # ftdata = []
        # for i in range(len(lonlats)):
            # ilon = lonlats.lon.values[i]; ilat =lonlats.lat.values[i];
            # idf = ydata[(ydata.lon==ilon)&(ydata.lat==ilat)]
            # its = fdata[(fdata.lon>=(ilon-0.25))&(fdata.lon<=(ilon+0.25))&(fdata.lat>=(ilat-0.25))&(fdata.lat<=(ilat+0.25))]
            # if (len(its)>0):
                # its = its.assign(date=lambda dfa: dfa.apply(lambda r: pd.date_range(r["disc_date"], r["end_date"]), axis=1)).explode("date")
                # zs = its.groupby(['date'],as_index=False).agg(
                        # {'fire_size': ['sum'],
                        # }) 
                # zs.columns = ['time','fire_size_FOD']
                # zs['lon'] =ilon; zs['lat'] =ilat; 
                # idf =idf.merge(zs,on =['lon', 'lat','time'],how='left')
            # else:
                # idf['fire_size_FOD']=np.nan
            # ##
            # ftdata +=[idf]
        # ##
        # ftdata=pd.concat(ftdata)    
        # ftdata.to_csv(os.path.join(indir,'HW_Wildfire_Drought_data',outfile),index=False)
    
# for y in years:
    # #ydata = tdf[tdf.time.dt.year==y]
    # #fdata = fod[fod.fire_year==y]
    # #fdata['discovery_date']=pd.to_datetime(fdata['discovery_date']);
    # outfile = 'CPC_Global_Unified_Temperature_with_HW_labels_Fires_'+str(y)+'.csv'
    # ydata = pd.read_csv(os.path.join(indir,'HW_Wildfire_Drought_data','CPC_Global_Unified_Temperature_with_HW_labels_Fires_'+str(y)+'.csv'))#.iloc[:, 1:]
    # ydata= ydata.drop(['fire_size_FOD'],axis=1)
    # ydata.to_csv(os.path.join(indir,'HW_Wildfire_Drought_data','CPC_Global_Unified_Temperature_with_HW_labels_Fires_'+str(y)+'.csv'),index=False)
    
    # nas = tdf.groupby(['lon','lat'])['tmax','tmin','tmean'].apply(lambda x: x.isnull().mean()).reset_index()
    # nas.to_csv(os.path.join(indir,'CPC_Global_Unified_Temperature_missing_value_summary.csv')) ## 24 grids with missing data
    #tdf=tdf.sort_values([])
    # hindex = ['HI02','HI04','HI05','HI06','HI09','HI10']
    # for y in years:
        # ydata = tdf[tdf.time.dt.year==y]
        # #ydata= tdf
        # lonlats = ydata.groupby(['lon', 'lat']).tmax.max().reset_index()
        # ldata = []
        # for i in range(len(lonlats)):
            # lon = lonlats.lon.values[i]; lat =lonlats.lat.values[i];
            # idata= ydata[(ydata.lon==lon)&(ydata.lat==lat)]
            # for h in range(len(hindex)):
                # hdata = pd.read_csv(os.path.join(indir,'Heatwave','Heatwave_'+hindex[h]+ '_conus.csv'))
                # hdata['start_date']=pd.to_datetime(hdata['start_date']);hdata['end_date']=pd.to_datetime(hdata['end_date'])
                # hdata = hdata[(hdata.lon==lon)&(hdata.lat==lat)]
                # hdata=hdata[hdata.start_date.dt.year==y]
                # if len(hdata)>0:
                    # cols = ['start_date', 'end_date','duration']
                    # df_merge = idata.merge(hdata[cols].drop_duplicates(),how='cross')   
                    # df_out = df_merge[df_merge['time'].between(df_merge['start_date'], df_merge['end_date'])]
                    # df_out[hindex[h]] =1
                    # #df_out=df_out[['lat', 'lon', 'time','start_date','end_date', 'duration', hindex[h]]]
                    # #df_out.columns = ['lat', 'lon', 'time',hindex[h]+'_start_date',hindex[h]+'_end_date', hindex[h]+'_duration', hindex[h]]
                    # df_out=df_out[['lat', 'lon', 'time', 'duration', hindex[h]]]
                    # df_out.columns = ['lat', 'lon', 'time', hindex[h]+'_duration', hindex[h]]
                    # idata=idata.merge(df_out, on = ['lat', 'lon', 'time'],how='left')
                    # idata[hindex[h]][idata[hindex[h]].isna()]=0;
                    # idata[hindex[h]+'_duration'][idata[hindex[h]+'_duration'].isna()]=0
                # else:
                    # idata[hindex[h]]=0;idata[hindex[h]+'_duration']=0;
            # ##
            # ldata += [idata]
        # ##
        # ldata= pd.concat(ldata)
        # ldata.to_csv(os.path.join(indir,'HW_Wildfire_Flood_data','CPC_Global_Unified_Temperature_with_HW_labels_'+str(y)+'.csv'))
            
 
    
# hdata = pd.read_csv(os.path.join(indir,'Heatwave','Heatwave_HI07_PNW.csv'))
# hdata['duration']=hdata.duration/np.timedelta64(1, 'D')
# hdata['start_date'] = pd.to_datetime(hdata['start_date'], format='%Y-%m-%d') #format='%Y-%m-%d' '%m/%d/%Y'
# hdata['year'] =hdata['start_date'].dt.year; hdata['year']=hdata['year'].astype('int')
# hdata['month'] =hdata['start_date'].dt.month
# h07f = hdata.groupby(['state','year'])['duration'].sum().reset_index()
# h07a = hdata.groupby(['state','year'])['tmax'].max().reset_index()
# fig, axes = plt.subplots(nrows=2, ncols=1,figsize=(6, 6))
# g1 =sns.lineplot(x='year',y='duration',hue='state',data=h07f,ax=axes[0])
# #g1.set_xticks(range(len(years))) # <--- set the ticks first
# #g1.set_xticklabels(labels=years)
# g1.set(xticks=np.arange(2001,2020,2))
# axes[0].legend(loc='upper left',ncol=3, title=None)
# axes[0].set_xlabel('Year')
# axes[0].set_ylabel('Total Heatwave days')
# g2=sns.lineplot(x='year',y='tmax',hue='state',data=h07a,ax=axes[1])
# #g2.set_xticks(range(len(years))) # <--- set the ticks first
# #g2.set_xticklabels(labels=years)
# g2.set(xticks=np.arange(2001,2020,2))
# axes[1].legend(loc='upper left',ncol=3, title=None)
# axes[1].set_xlabel('Year')
# axes[1].set_ylabel('Highest temperature (c)')
# plt.tight_layout()
# plt.savefig(os.path.join('/home/linx882/LDRD_extremes/Heatwave/','HW'+'_HI07_tsplot.png'), bbox_inches='tight',dpi=300)


# hdata = pd.read_csv(os.path.join(indir,'Heatwave','Heatwave_HI10_PNW.csv'))
# hdata['duration']=hdata.duration/np.timedelta64(1, 'D')
# hdata['start_date'] = pd.to_datetime(hdata['start_date'], format='%Y-%m-%d') #format='%Y-%m-%d' '%m/%d/%Y'
# hdata['year'] =hdata['start_date'].dt.year; hdata['year']=hdata['year'].astype('int')
# hdata['month'] =hdata['start_date'].dt.month
# h10f = hdata.groupby(['state','year'])['duration'].sum().reset_index()
# h10a = hdata.groupby(['state','year'])['tmin'].max().reset_index()
# fig, axes = plt.subplots(nrows=2, ncols=1,figsize=(6, 6))
# g1 =sns.lineplot(x='year',y='duration',hue='state',data=h10f,ax=axes[0])
# #g1.set_xticks(range(len(years))) # <--- set the ticks first
# #g1.set_xticklabels(labels=years)
# g1.set(xticks=np.arange(2001,2020,2))
# axes[0].legend(loc='upper left',ncol=3, title=None)
# axes[0].set_xlabel('Year')
# axes[0].set_ylabel('Total Heatwave days')
# g2=sns.lineplot(x='year',y='tmin',hue='state',data=h10a,ax=axes[1])
# #g2.set_xticks(range(len(years))) # <--- set the ticks first
# #g2.set_xticklabels(labels=years)
# g2.set(xticks=np.arange(2001,2020,2))
# axes[1].legend(loc='upper left',ncol=3, title=None)
# axes[1].set_xlabel('Year')
# axes[1].set_ylabel('Highest Tmin (c)')
# plt.tight_layout()
# plt.savefig(os.path.join('/home/linx882/LDRD_extremes/Heatwave/','HW'+'_HI10_tsplot.png'), bbox_inches='tight',dpi=300)


# hdata = pd.read_csv(os.path.join(indir,'Heatwave','Heatwave_HI06_PNW.csv'))
# hdata['duration']=hdata.duration/np.timedelta64(1, 'D')
# hdata['start_date'] = pd.to_datetime(hdata['start_date'], format='%Y-%m-%d') #format='%Y-%m-%d' '%m/%d/%Y'
# hdata['year'] =hdata['start_date'].dt.year; hdata['year']=hdata['year'].astype('int')
# hdata['month'] =hdata['start_date'].dt.month
# hdata['ym'] = hdata['start_date'].apply(lambda x: x.strftime('%Y-%m')) 
# #idf = hdata[(hdata.lon==-121.75)&(hdata.lat==44.75)]
# idf = hdata.groupby(['lat', 'lon','state','year']).duration.sum().reset_index()
# #idf = idf.groupby(['year']).tmax.max().reset_index()



# # v=v.drop('state',axis=1)
# # v = v.merge(cp[['cp','state']],on='cp',how='left')
# # v= v.dropna(subset=['state'],axis=0)

# PNW= ['Washington','Idaho','Oregon']
# sv3 = vp[vp.state.isin(PNW)]
# sv3['dtime'] = pd.to_datetime(sv3['dtime'], format='%Y-%m-%d') #format='%Y-%m-%d' '%m/%d/%Y'
# sv3['year'] =sv3['dtime'].dt.year; sv3['year']=sv3['year'].astype('int')
# sv3['month'] =sv3['dtime'].dt.month; 
# sv3['ym'] = sv3['dtime'].apply(lambda x: x.strftime('%Y-%m')) 
# sgls = sv3.groupby(['lat', 'lon']).gl.count().reset_index()

# yts = sv3.groupby(['lat', 'lon','year']).gl.count().reset_index()
# # yts = sv3.groupby(['lat', 'lon','year']).agg(
            # # {'gl': ['count'],
             # # 'tmax': ['max'],
            # # })
# #yts[(yts.gl>5)&(yts.gl<10)]
# data = yts.copy()
# # its = data[(data.lon>=-122)&(data.lon<=-121.5)&(data.lat>=44.5)&(data.lat<=45)]
# # its = its.groupby(['year']).gl.sum().reset_index()
# # its =its.merge(idf,on='year',how='outer')

# #idf = hdata.groupby(['lat', 'lon','state','year']).duration.sum().reset_index()
# df = hdata.groupby(['lat', 'lon','state','year']).agg(
            # {'duration': ['sum'],
             # 'tmax': ['max'],
            # }).reset_index()
# df.columns = ['lat', 'lon','state','year', 'duration','tmax']
# mdata =[]
# for i in range(len(df)):
    # ilon = df.lon[i]; ilat = df.lat[i]; state =df.state[i]
    # idf = df[(df.lon==ilon)&(df.lat==ilat)]
    # its = data[(data.lon>=(ilon-0.25))&(data.lon<=(ilon+0.25))&(data.lat>=(ilat-0.25))&(data.lat<=(ilat+0.25))]
    # its = its.groupby(['year']).gl.sum().reset_index()
    # its['lon'] =ilon; its['lat'] =ilat; its['state'] =state;
    # idf =idf.merge(its,on =['lon', 'lat','year','state'],how='outer')
    # mdata += [idf.copy()]
    
# mdata=pd.concat(mdata)
# #mdata.to_csv(os.path.join(indir,'HW_fires_annual.csv'))
# mdata = pd.read_csv(os.path.join(indir,'HW_fires_annual.csv'))
# sdata = mdata.drop_duplicates()
# idata = sdata[(sdata.lon==-121.75)&(sdata.lat==44.75)]
# idata=idata.sort_values(['year'])
# cc = sdata.groupby(['lat', 'lon'])[['duration','gl']].corr().unstack().iloc[:,1].reset_index()
# cc.columns = ['lat', 'lon','corr']
# cc[cc['corr']>0.7]

# data=cc
# fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(6, 6))
# #ax = conus_new.plot(facecolor="none",edgecolor=conus['region'],linewidth=2)
# ax.set_title('Correlation between fire grids and HW days',fontsize=8)
# divider = make_axes_locatable(ax)
# cax = divider.append_axes("right", size="3%", pad=0.01)
# x, y = data['lon'].values, data['lat'].values
# # x, y = tdata_ngp['lon'].values, tdata_ngp['lat'].values
# im=ax.scatter(-123.75, 46.75, c='#EE3B3B', marker="*",s=25,alpha=0.6) #data[(data.lon==-122.75)&(data.lat==47.75)]['corr'],marker="s",alpha=0.8, cmap=color_ramp
# im=ax.scatter(-120.75, 45.25, c='#EE3B3B', marker="*",s=25,alpha=0.6) #data[(data.lon==-117.75)&(data.lat==44.75)]['corr']
# im=ax.scatter(-115.25, 46.25, c='#EE3B3B', marker="*",s=25,alpha=0.6) #data[(data.lon==-111.75)&(data.lat==42.25)]['corr']
# im=ax.scatter(x, y, c=data['corr'], marker="s",s=15,cmap=plt.get_cmap("coolwarm"),alpha=0.6) #,marker="s",alpha=0.8, cmap=color_ramp
# #im=ax.scatter(x, y, marker="*",s=0.01) #,marker="s"
# #im.set_clim(data['gl'].quantile(.1),data['gl'].quantile(.9))
# #im.set_clim(1,3)
# #conus.plot(facecolor="none",edgecolor='r',linewidth=1,alpha=0.2,ax=ax)
# PNW_shape.plot(facecolor="none",edgecolor='k',linewidth=1,ax=ax)
# fig.colorbar(im, cax=cax, orientation='vertical',shrink=0.9)
# ax.set_xlabel('Lon')
# ax.set_ylabel('Lat')
# plt.savefig(os.path.join('/home/linx882/LDRD_extremes/','fire_plots','annual_grid_duration_s3.png'),
            # transparent = False,facecolor = 'white',bbox_inches='tight',dpi=300)
            

# kdata =sdata[(sdata.lon==-123.75)&(sdata.lat==46.75)]
# kdata = sdata[(sdata.lon==-120.75)&(sdata.lat==45.25)]
# kdata = sdata[(sdata.lon==-115.25)&(sdata.lat==46.25)]
# #
# fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(4, 4))
# ax=sns.regplot(x=kdata.duration, y=kdata.gl,label=True, ci=None)
# # label points on the plot
# for x, y,z in zip(kdata['duration'], kdata['gl'],kdata['year']):
    # # the position of the data label relative to the data point can be adjusted by adding/subtracting a value from the x &/ y coordinates
    # ax.text(x = x, # x-coordinate position of data label
    # y = y-0.2, # y-coordinate position of data label, adjusted to be 150 below the data point
    # s = '{:.0f}'.format(z), # data label, formatted to ignore decimals
    # color = 'purple') # set colour of line
# ax.set_xlabel('HW days')
# ax.set_ylabel('Annual total grids')
# plt.savefig(os.path.join('/home/linx882/LDRD_extremes/','fire_plots','scatter_grid_duration_WA.png'),
            # transparent = False,facecolor = 'white',bbox_inches='tight',dpi=300)
# plt.show()




# cc = sdata.groupby(['lat', 'lon'])[['tmax','gl']].corr().unstack().iloc[:,1].reset_index()
# cc.columns = ['lat', 'lon','corr']

# data=cc
# fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(6, 6))
# #ax = conus_new.plot(facecolor="none",edgecolor=conus['region'],linewidth=2)
# ax.set_title('Correlation between fire grids and Max temperature',fontsize=8)
# divider = make_axes_locatable(ax)
# cax = divider.append_axes("right", size="3%", pad=0.01)
# x, y = data['lon'].values, data['lat'].values
# # x, y = tdata_ngp['lon'].values, tdata_ngp['lat'].values
# im=ax.scatter(x, y, c=data['corr'], marker="s",s=15,cmap=plt.get_cmap("coolwarm"),alpha=0.6) #,marker="s",alpha=0.8, cmap=color_ramp
# #im=ax.scatter(x, y, marker="*",s=0.01) #,marker="s"
# #im.set_clim(data['gl'].quantile(.1),data['gl'].quantile(.9))
# #im.set_clim(1,3)
# #conus.plot(facecolor="none",edgecolor='r',linewidth=1,alpha=0.2,ax=ax)
# PNW_shape.plot(facecolor="none",edgecolor='k',linewidth=1,ax=ax)
# fig.colorbar(im, cax=cax, orientation='vertical',shrink=0.9)
# ax.set_xlabel('Lon')
# ax.set_ylabel('Lat')
# plt.savefig(os.path.join('/home/linx882/LDRD_extremes/','fire_plots','annual_grid_tmax.png'),
            # transparent = False,facecolor = 'white',bbox_inches='tight',dpi=300)
# ########################################
# ## get percentiles for each grid
# ftmins = []
# cdata =[]
# for y in years:
    # data = xr.open_dataset(os.path.join(tdir,[f for f in ftmins if str(y) in f][0]))
    # # tmin.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
    # # tmin.rio.write_crs("EPSG:4326", inplace=True) #data.rio.set_crs("EPSG:4326") 
    # #tmin = xr.open_dataset(os.path.join(tdir,[f for f in ftmins if str(year) in f][0]))
    # # tmin.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
    # # tmin.rio.write_crs("EPSG:4326", inplace=True) #data.rio.set_crs("EPSG:4326") 
    # data = data.assign_coords(lon=(((data.lon + 180) % 360) - 180))
    # #data.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
    # sdata = data.where((data.lon>=-125)&(data.lon<=-66)&(data.lat>=25)&(data.lat<=50), drop=True) #-125,25,-66,50
    # sdata.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
    # sdata.rio.write_crs("EPSG:4326", inplace=True)
    # sdata=  sdata.rio.clip(conus.geometry.apply(mapping), conus.crs, drop=True)
    # cdata +=[sdata]

# kcdata= xr.concat(cdata,dim ='time')
# qt_dims = ("time")
# qt_values = (0.75,0.81, 0.9,0.95,0.975, 0.98, 0.99)
# isdf = kcdata.quantile(qt_values, dim=qt_dims)
# #tmin['Date'] = str(y)+'-'+tmin.Month.apply(str).str.zfill(2)+'-'+tmin.Day.apply(str).str.zfill(2)
# isdf = isdf.to_dataframe().reset_index(); isdf = isdf.dropna(); 
# #isdf=isdf.drop('spatial_ref',axis=1)
# #isdf['quantile'] = isdf['quantile']*100#.astype(int)
# data = isdf.pivot(index=['lon','lat'], columns='quantile', values='tmin')
# data.columns = ['Tmin75', 'Tmin81', 'Tmin90','Tmin95', 'Tmin97.5', 'Tmin98', 'Tmin99']
# data= data.reset_index()
# data.to_csv(os.path.join(tdir,'Tmin_thresholds.csv'))

# isdf['state'] =state

# pmax = pd.read_csv(os.path.join(tdir,'Tmax_thresholds.csv'))
# pmin = pd.read_csv(os.path.join(tdir,'Tmin_thresholds.csv'))
# pdata = pmax.merge(pmin,on=['lon', 'lat'])
# sdata =[]
# for s in states:
    # idata = point_in_state(shapes,s,pdata)
    # sdata+=[idata]
# sdata=pd.concat(sdata)
# sdata.to_csv(os.path.join(tdir,'Tmax_Tmin_thresholds.csv'),index=False)






# data=isdf
# fig, ax = plt.subplots(figsize=(6,4),dpi=300)
# #ax = conus_new.plot(facecolor="none",edgecolor=conus['region'],linewidth=2)
# ax.set_title('Fire Grids',fontsize=12)
# x, y = data['lon'].values, data['lat'].values
# # x, y = tdata_ngp['lon'].values, tdata_ngp['lat'].values
# im=ax.scatter(x, y, c=data.tmin, marker="s",s=2,cmap=plt.get_cmap("jet"),alpha=0.2) #,marker="s"
# #im=ax.scatter(x, y, marker="*",s=0.01) #,marker="s"
# #im.set_clim(-0.05,0.05)
# #conus.plot(facecolor="none",edgecolor='r',linewidth=1,alpha=0.2,ax=ax)
# shape.plot(facecolor="none",edgecolor='r',linewidth=1,alpha=0.2,ax=ax)
# cbar = plt.colorbar(im, shrink=0.65, aspect=20, fraction=.12,pad=.02)
# #ax = conus_shapes.plot(facecolor="none",edgecolor='k',linewidth=3)
# # cbar.set_label('Maximum 3-h accumulation',size=10)
# # access to cbar tick labels:
# cbar.ax.tick_params(labelsize=10) 
# plt.ylabel('Lat',fontsize=12)
# plt.xlabel('Lon',fontsize=12)
# plt.xticks(fontsize=12)
# plt.yticks(fontsize=12)
# plt.show()