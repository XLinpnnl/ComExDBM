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


# def GUTemp_state_daily_tmax_tmin(tdir,year,state,shapes):
    # # read in CPC Global Unified Temperature(e.g.,tmin,tmax) by state
    # # calculate daily men (tmean = (tmax+tmin)/2)
    # ftmins = [fn for fn in os.listdir(tdir) if 'tmin' in fn]
    # ftmins.sort()
    # ftmaxs = [fn for fn in os.listdir(tdir) if 'tmax' in fn]
    # ftmaxs.sort()
    # shape = shapes[shapes.NAME==state]
    # bounds = shape.geometry.bounds
    # tmax = xr.open_dataset(os.path.join(tdir,[f for f in ftmaxs if str(year) in f][0]))
    # # tmin.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
    # # tmin.rio.write_crs("EPSG:4326", inplace=True) #data.rio.set_crs("EPSG:4326") 
    # tmin = xr.open_dataset(os.path.join(tdir,[f for f in ftmins if str(year) in f][0]))
    # # tmin.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
    # # tmin.rio.write_crs("EPSG:4326", inplace=True) #data.rio.set_crs("EPSG:4326") 
    # data = xr.merge([tmax,tmin])
    # data['tmean'] = (data['tmax']+data['tmin'])/2
    # data = data.assign_coords(lon=(((data.lon + 180) % 360) - 180))
    # #data.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
    # sdata = data.where((data.lon>=bounds.minx.values)&(data.lon<=bounds.maxx.values)&(data.lat>=bounds.miny.values)&(data.lat<=bounds.maxy.values), drop=True) #-125,25,-66,50
    # sdata.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
    # sdata.rio.write_crs("EPSG:4326", inplace=True)
    # sdata0=  sdata.rio.clip(shape.geometry.apply(mapping), shape.crs, drop=True)
    # #tmin['Date'] = str(y)+'-'+tmin.Month.apply(str).str.zfill(2)+'-'+tmin.Day.apply(str).str.zfill(2)
    # isdf = sdata.to_dataframe().reset_index();    
    # isdf = isdf.dropna(); isdf=isdf.drop('spatial_ref',axis=1)
    # isdf = isdf.rename(columns={'time': 'date'})
    # isdf['state'] =state
    # return isdf

def GUTemp_state_daily_temp(data,state,shapes):
    # read in CPC Global Unified Temperature(e.g.,tmin,tmax) by state
    # calculate daily men (tmean = (tmax+tmin)/2)
    shape = shapes[shapes.NAME==state]
    bounds = shape.geometry.bounds
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


def get_HIxxt1(tdata,state,shapes,pdata,var,threshold=99,days=3,metric='HI04'):
    ydata = GUTemp_state_daily_temp(tdata,state,shapes)
    lonlats = ydata.groupby(['lon', 'lat']).tmax.max().reset_index()
    hdata = []
    for i in range(len(lonlats)):
        lon = lonlats.lon.values[i]; lat =lonlats.lat.values[i];
        Tm= pdata[var+str(threshold)][(pdata.lon==lon)&(pdata.lat==lat)].values[0]  ## T1>90 percentiles #T2>75 percentiles
        #Tm =Tm*9/5+32;
        sdata = ydata[(ydata.lon==lon)&(ydata.lat==lat)];sdata= sdata.sort_values(by=['date'])#sdata = ydata[['Date', 'FIPS', 'T2_mean', 'T2_max', 'T2_min']][ydata.FIPS==ci]
        #sdata.date = sdata['date'].apply(pd.to_datetime)
        sdata['HT'] = [1 if t >Tm else 0 for t in sdata[var]]
        #sdata['Hn']=sdata['HT'].groupby((sdata['HT'] != sdata['HT'] .shift()).cumsum()).cumsum()
        sdf = sdata[sdata.HT==1]
        sdf['grp_date'] = sdf.date.diff().dt.days.ne(1).cumsum()
        #sdf['state'] = state
        cols = ['tmax', 'tmin', 'tmean']
        zs = sdf.groupby(['lat', 'lon','grp_date'],as_index=False).agg(
                {cols[0]: ['max','mean'],
                 cols[1]: ['max','mean'],
                 cols[2]: ['max','mean'],
                })   
        columns=["_".join(x) for x in zs.columns.ravel()]
        columns = ['lat', 'lon','grp_date'] +columns[3:]; columns = [c.replace('_mean','_avg') for c in columns]
        zs.columns = columns
        z3 = sdf.reset_index().groupby(['lat', 'lon','grp_date'])['date'].agg([min, max])
        z3= z3.assign(duration=(z3['max'] - z3['min'])+ timedelta(days=1)).reset_index()
        #z3 = z3.reset_index()
        #result=z1.merge(z2, on=['FIPS','grp_date'])
        result=reduce(lambda x, y: pd.merge(x, y, on = ['lat', 'lon','grp_date']), [zs,z3])
        result = result.rename(columns={'min': 'start_date','max': 'end_date'})
        result[var+str(threshold)] = Tm
        hdata += [result]
    #
    ##
    hdata= pd.concat(hdata, ignore_index=True)
    hdata=hdata[hdata.duration>=timedelta(days=days)]
    #hdata.to_csv(os.path.join('/home/linx882/TGW/','HW_data_'+senario,metric+"_Jul_Aug"+str(year)+".csv"),index=False)
    return hdata


def get_HIxxt2(tdata,state,shapes,pdata,var,t1=90,t2=75,metric='HI07'):
    #pdata.FIPS=pdata.FIPS.astype(str).str.zfill(5)
    #cid = ydata.FIPS.unique()
    #cid = pdata.FIPS.unique()
    ydata = GUTemp_state_daily_temp(tdata,state,shapes)
    lonlats = ydata.groupby(['lon', 'lat']).tmax.max().reset_index()
    #ydata =pdata[pdata.state==state]
    hdata = []
    for i in range(len(lonlats)):
        #(1°C × 9/5) + 32 = 33.8°F
        lon = lonlats.lon.values[i]; lat =lonlats.lat.values[i];
        TT= pdata[(pdata.lon==lon)&(pdata.lat==lat)]
        #T1= TT[var.capitalize()+str(t1)].values[0] ; T2= TT[var.capitalize()+str(t2)].values[0] ## T1>90 percentiles #T2>75 percentiles
        T1= TT[var+str(t1)].values[0] ; T2= TT[var+str(t2)].values[0] ## T1>90 percentiles #T2>75 percentiles
        #T1 =T1*9/5+32; T2 = T2*9/5+32;
        sdata = ydata[(ydata.lon==lon)&(ydata.lat==lat)];sdata= sdata.sort_values(by=['date'])
        if len(sdata)>0:
            #sdata.Date = sdata['Date'].apply(pd.to_datetime)
            sdata['T1'] = [1 if t >T1 else 0 for t in sdata[var]]
            sdata['T2'] = [1 if t >T2 else 0 for t in sdata[var]]
            #sdata['Hn']=sdata['HT'].groupby((sdata['HT'] != sdata['HT'] .shift()).cumsum()).cumsum()
            #sdf = sdata[(sdata.T1==1)]
            sdf = sdata[(sdata.T2==1)]
            if len(sdf)>0:
                sdf['grp_date'] = sdf.date.diff().dt.days.ne(1).cumsum()
                sdf=sdf.groupby('grp_date').filter(lambda x : len(x)>=3)
                HT2=sdf['grp_date'].unique()
                idata =[]
                for ht in HT2:
                    isdf= sdf[sdf.grp_date==ht]
                    kdf= isdf[(isdf.T1==1)]
                    if len(kdf)>=3:
                        kdf['grp_date2'] = kdf.date.diff().dt.days.ne(1).cumsum()
                        kdf=kdf.groupby('grp_date2').filter(lambda x : len(x)>=3)
                        #td =kdf.groupby(['grp_date2']).grp_date2.count().max()
                        if (len(kdf)>=3):
                            tm =isdf['tmax'].mean()
                            while tm<T1:
                                isdf=isdf.drop(isdf[var].idxmin()) #isdf[-isdf.T2_max==isdf.T2_max.min()]
                                isdf['grp_date'] = isdf.date.diff().dt.days.ne(1).cumsum()
                                isdf=isdf.groupby('grp_date').filter(lambda x : len(x)>=3)
                                itd =isdf.groupby(['grp_date'])[var].mean().reset_index()
                                itd1 =itd[itd.tmax>=T1]
                                if (len(itd1)>0):
                                    isdf =isdf[isdf.grp_date.isin(itd1.grp_date)]
                                    tm = itd1.tmax.min()
                                else:
                                    tm = itd.tmax.max()
                            # #ig = itd.grp_date[itd[var]==itd[var].max()].values[0]
                            # isdf = isdf[isdf.grp_date==ig]
                            # tm =isdf[var].mean()
                            else :
                                isdf=isdf
                            ##
                            isdf['grp_date'] = str(ht)+'_'+isdf['grp_date'].apply(str)
                            #isdf['NAME'] = pdata.NAME[pdata.FIPS==ci].values[0]
                            #z1 = isdf.groupby(['lat', 'lon','state','grp_date'])['tmax'].max().reset_index()
                            #z2 = isdf.groupby(['lat', 'lon','state','grp_date'])['tmax'].mean().reset_index();z2.columns = ['lat', 'lon','state','grp_date',var+'_avg']
                            cols = ['tmax', 'tmin', 'tmean']
                            zs = isdf.groupby(['lat', 'lon','grp_date'],as_index=False).agg(
                                    {cols[0]: ['max','mean'],
                                     cols[1]: ['max','mean'],
                                     cols[2]: ['max','mean'],
                                    })   
                            columns=["_".join(x) for x in zs.columns.ravel()]
                            columns = ['lat', 'lon','grp_date'] +columns[3:]; columns = [c.replace('_mean','_avg') for c in columns]
                            zs.columns = columns
                            z3 = isdf.reset_index().groupby(['lat', 'lon','grp_date'])['date'].agg([min, max])
                            z3 = z3.assign(duration=(z3['max'] - z3['min'])+ timedelta(days=1)).reset_index()
                            #z3 = z3.reset_index()
                            #result=z1.merge(z2, on=['FIPS','grp_date'])
                            iresult=reduce(lambda x, y: pd.merge(x, y, on = ['lat', 'lon','grp_date']), [zs,z3])
                            iresult = iresult.rename(columns={'min': 'start_date','max': 'end_date'})
                            iresult['T'+str(t1)] = T1; iresult['T'+str(t2)] = T2; 
                            idata += [iresult]
                #ht
                if len(idata)>0:
                    result = pd.concat(idata, ignore_index=True)
                    hdata += [result]
    ##
    if len(hdata)>0:
        hdata= pd.concat(hdata, ignore_index=True)
        #hdata=hdata[hdata.duration>=timedelta(days=days)]
        #hdata.to_csv(os.path.join('/home/linx882/TGW/','HW_data_'+senario,metric+"_Jul_Aug"+str(year)+".csv"),index=False)
        return hdata

def plot_HW(hdata,var,cshape, note):
    if hdata.duration.dtype!=np.float64:
        hdata['duration']=hdata.duration/np.timedelta64(1, 'D')
    ##
    hws = hdata[hdata.duration>=3]
    #The sum of HW days (HWF)
    hwf = hws.groupby(['lat', 'lon'])['duration'].sum().reset_index()
    hwf.columns =['lat', 'lon', 'HWF']
    # The highest temperature of the hottest event (HWA) 
    hwa = hws.groupby(['lat', 'lon'])['tmax_max'].max().reset_index()
    hwa.columns =['lat', 'lon', 'HWA']
    #The spatial extent of HW. 
    #spe = len(hws.FIPS.unique())/len(hdata.FIPS.unique())
    data = reduce(lambda left,right: pd.merge(left,right,on=['lat', 'lon'],
                                                how='left'), [hwf,hwa])
    #hwshapes=hwshapes.fillna(0)
    ramp_colors=['#3C4CA7','#486ABB','#77D2E3','#5CC58A','#3CB73D','#5CC133','#B7DC17','#F6ED0B','#FABE09','#F47D13','#EE371F','#EC1C24','#A0164B','#4F2783']
    color_ramp = LinearSegmentedColormap.from_list( 'my_list', [ Color( c1 ).rgb for c1 in ramp_colors ] )
    ##
    fig, axes = plt.subplots(nrows=2, ncols=1,figsize=(6, 4))
    #ax = conus_new.plot(facecolor="none",edgecolor=conus['region'],linewidth=2)
    axes[0].set_title('Sum of '+note+' Positive Days',fontsize=8)
    divider = make_axes_locatable(axes[0])
    cax = divider.append_axes("right", size="3%", pad=0.01)
    x, y = data['lon'].values, data['lat'].values
    # x, y = tdata_ngp['lon'].values, tdata_ngp['lat'].values
    im=axes[0].scatter(x, y, c=data['HWF'], marker="s",s=10,cmap=plt.get_cmap("jet"),alpha=0.8) #,marker="s",alpha=0.8, cmap=color_ramp
    #im=ax.scatter(x, y, marker="*",s=0.01) #,marker="s"
    im.set_clim(data['HWF'].quantile(.1),data['HWF'].quantile(.9))
    #conus.plot(facecolor="none",edgecolor='r',linewidth=1,alpha=0.2,ax=ax)
    cshape.plot(facecolor="none",edgecolor='k',linewidth=1,ax=axes[0])
    fig.colorbar(im, cax=cax, orientation='vertical',shrink=0.9)
    axes[0].set_xlabel('Lon')
    axes[0].set_ylabel('Lat')
    #cbar = plt.colorbar(im, shrink=0.9, aspect=20, fraction=.12,pad=.02)
    #ax = conus_shapes.plot(facecolor="none",edgecolor='k',linewidth=3)
    # cbar.set_label('Maximum 3-h accumulation',size=10)
    # access to cbar tick labels:
    #cbar.ax.tick_params(labelsize=10) 
    axes[1].set_title('Highest '+var+'(C)',fontsize=8)
    divider = make_axes_locatable(axes[1])
    cax = divider.append_axes("right", size="3%", pad=0.01)
    x, y = data['lon'].values, data['lat'].values
    # x, y = tdata_ngp['lon'].values, tdata_ngp['lat'].values
    im2=axes[1].scatter(x, y, c=data['HWA'], marker="s",s=10,cmap=plt.get_cmap("jet"),alpha=0.8) #,marker="s",alpha=0.8, cmap=color_ramp
    #im=ax.scatter(x, y, marker="*",s=0.01) #,marker="s"
    im2.set_clim(data['HWA'].quantile(.1),data['HWA'].quantile(.9))
    #conus.plot(facecolor="none",edgecolor='r',linewidth=1,alpha=0.2,ax=ax)
    cshape.plot(facecolor="none",edgecolor='k',linewidth=1,ax=axes[1])
    fig.colorbar(im2, cax=cax, orientation='vertical',shrink=0.9)
    #cbar2 = plt.colorbar(im2, shrink=0.9, aspect=20, fraction=.12,pad=.02)
    #ax = conus_shapes.plot(facecolor="none",edgecolor='k',linewidth=3)
    # cbar.set_label('Maximum 3-h accumulation',size=10)
    # access to cbar tick labels:
    #cbar2.ax.tick_params(labelsize=10) 
    axes[1].set_xlabel('Lon')
    axes[1].set_ylabel('Lat')
    plt.tight_layout()
    plt.pause(1e-13)
    #plt.subplots_adjust(left=0.125,bottom=0.05, right=0.95, top=0.95, wspace=0.15, hspace=0.25)
    plt.savefig(os.path.join('/home/linx882/LDRD_extremes/Heatwave/','HW'+'_'+note+'_metrics.png'), bbox_inches='tight',dpi=300)
    return hwf,hwa

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
    PNW_shape = conus[conus.NAME.isin(['Washington','Idaho','Oregon'])]
    states = conus.NAME.unique().tolist()
    PNW= ['Washington','Idaho','Oregon']
    #pdata =pd.read_csv(os.path.join(tdir,'Tmax_Tmin_thresholds.csv'))
    pdata= pd.read_csv(os.path.join(tdir,'Temp_thresholds.csv'))
    #
    # hdata =[]
    # for y in years:
         # for s in PNW:
            # idata = get_HIxxt2(year=y,state=s,shapes=conus,pdata=pdata,var='tmax',t1=90,t2=75,metric='HI07')
            # hdata += [idata]
    # #
    # hdata=pd.concat(hdata)
    # hdata.to_csv(os.path.join(indir,'Heatwave','Heatwave_HI07_PNW.csv'),index=False)
    
    # hdata = pd.read_csv(os.path.join(indir,'Heatwave','Heatwave_HI07_PNW.csv'))
    # plot_HW(hdata=hdata,var='tmax',cshape=PNW_shape, note='HI07')
    # hdata =[]
    # for y in years:
         # for s in PNW:
            # idata = get_HIxxt2(year=y,state=s,shapes=conus,pdata=pdata,var='tmax',t1=97.5,t2=81,metric='HI06')
            # hdata += [idata]
    # #
    # hdata=pd.concat(hdata)
    # hdata.to_csv(os.path.join(indir,'Heatwave','Heatwave_HI06_PNW.csv'),index=False)
    # hdata = pd.read_csv(os.path.join(indir,'Heatwave','Heatwave_HI06_PNW.csv'))
    # plot_HW(hdata=hdata,var='tmax',cshape=PNW_shape, note='HI06')
    tdata = xr.open_dataset(os.path.join(tdir,'tmax_tmin_tmean_conus_2001_2020.nc'))
    hindex = ['HI02','HI04','HI05','HI09']
    tvars = ['tmean','tmean','tmax','tmin']
    thresholds = [95, 99,95,95]
    #
    for i in range(len(hindex)):
        hdata =[]
        for s in states:
            idata = get_HIxxt1(tdata=tdata,state=s,shapes=conus,pdata=pdata,var=tvars[i],threshold=thresholds[i],days=3,metric=hindex[i])
            hdata += [idata]
        #
        hdata=pd.concat(hdata)
        hdata.to_csv(os.path.join(indir,'Heatwave','Heatwave_'+hindex[i]+ '_conus.csv'),index=False)
    #
    hindex2 = ['HI06','HI10']
    tvars2 = ['tmax','tmin']
    #
    for i in range(len(hindex2)):
        hdata =[]
        for s in states:
            idata = get_HIxxt2(tdata=tdata,state=s,shapes=conus,pdata=pdata,var=tvars2[i],t1=97.5,t2=81,metric=hindex2[i])
            hdata += [idata]
        #
        hdata=pd.concat(hdata)
        hdata.to_csv(os.path.join(indir,'Heatwave','Heatwave_'+hindex2[i]+ '_conus.csv'),index=False)
    
    #
    hindex = ['HI02','HI04','HI05','HI06','HI09','HI10']
    tvars = ['tmean','tmean','tmax','tmax','tmin','tmin']
    for h in range(len(hindex)):
        hdata = pd.read_csv(os.path.join(indir,'Heatwave','Heatwave_'+hindex[h]+ '_conus.csv'))
        hdata=hdata.drop_duplicates()
        hdata.to_csv(os.path.join(indir,'Heatwave','Heatwave_'+hindex[h]+ '_conus.csv'),index=False)
        plot_HW(hdata=hdata,var=tvars[h],cshape=conus, note=hindex[h]+'_conus') 

    
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
# #################
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