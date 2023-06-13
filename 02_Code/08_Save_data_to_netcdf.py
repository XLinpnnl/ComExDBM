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
import itertools
import warnings
warnings.filterwarnings("ignore")
##############################################################################
indir = '/qfs/projects/comexdbm/LDRD_extremes/'
conus = gpd.read_file(os.path.join(indir,'CONUS4326','contig_us.shp'), crs="EPSG:4326")
years= list(range(2001, 2021))

ndir ='/qfs/projects/comexdbm/LDRD_extremes/NARR_vars_daily/NARR_regrid/'
nfiles = [f for f in os.listdir(ndir)if 'csv' in f]

mdir ='/qfs/projects/comexdbm/LDRD_extremes/MERRA2_vars_daily/MERRA_regrid/'
mfiles = [f for f in os.listdir(mdir)if 'csv' in f]

hindex = ['HI02','HI04','HI05','HI06','HI09','HI10']
findex= ['FHS_c9','FHS_c8c9','maxFRP_max_c9','BA_km2']
dindex1 = ['pdsi','spi14d','spi30d','spi90d','spei14d','spei30d','spei90d']
dindex2 =[d+'_category' for d in dindex1]
dindex = dindex1+dindex2

tdata = xr.open_dataset(os.path.join(indir,'CPC_Global_Unified_Temperature','tmax.2001.nc'))
pdata = xr.open_dataset(os.path.join(indir,'CPC_Global_Unified_Temperature','precip.V1.0.2001.nc'))

for y in years:
    #ydata = tdf[tdf.time.dt.year==y]
    yfile = 'CPC_Global_Unified_Temperature_HFD_'+str(y)+'.csv'
    #files= [fn for fn in os.listdir(os.path.join(indir,'HW_Wildfire_Drought_data'))]
    ydata = pd.read_csv(os.path.join(indir,'HW_Wildfire_Drought_data',yfile))#.iloc[:, 1:]
    #ydata=ydata.rename({'maxFRP_max_c9':'maxFRP_max'})
    hcols =['date', 'lon', 'lat','tmax', 'tmin', 'tmean']+[c for c in ydata.columns if c in hindex]
    fcols =['date', 'lon', 'lat']+[c for c in ydata.columns if c in findex]
    dcols =['date', 'lon', 'lat']+[c for c in ydata.columns if c in dindex]
    #
    hds = ydata[hcols].set_index(['date', 'lon', 'lat']).to_xarray()
    hds.attrs = {'title': 'CPC Global Temperature and Heatwave Labels', 
                'description': 'Gridded daily temperature(min/max/mean), and heatwave labels uding different definitions',
                'Conventions': 'CF-1.2','version': 'V1.0', 
                'SpatialCoverage': 'contiguous United States', 'SpatialResolution': '0.5 x 0.5', 
                'TemporalRange': str(y)+'-01-01 -> '+str(y)+'-12-31',  'TemporalResolution': 'daily', 
                'Format': 'NetCDF-4/HDF-5', 
                'comment': 'the original temperature data are originally created 9/2016 by CAS NOAA/ESRL PSD',
                'Original Tepmeattue Data Source': 'ftp://ftp.cpc.ncep.noaa.gov/precip/wd52ws/global_temp/', 
                'References': 'https://www.psl.noaa.gov/data/gridded/data.cpc.globaltemp.html'
                }
    hds.to_netcdf(os.path.join(indir,'ComExDBM/04_OutputData','01_HeatWave','ComExDBM_'+str(y)+'_HeatWave'+'_V01.nc'))
    #
    fds = ydata[fcols].set_index(['date', 'lon', 'lat']).to_xarray()
    fds.attrs = {'title': 'Fire related features from MODIS products', 'version': 'V1.0', 
                'description': 'Gridded daily fire data,including number of fire hot spot(i.e.,FHS_c9,FHS_c8c9), intensity(i.e.,maxFRP_max_c9 ) and bured areas(i.e.,BA_km2)',
                'SpatialCoverage': 'contiguous United States', 'SpatialResolution': '0.5 x 0.5', 
                'TemporalRange': str(y)+'-01-01 -> '+str(y)+'-12-31',  'TemporalResolution': 'daily', 
                'Format': 'NetCDF-4/HDF-5', 
                'comment': 'Fire features here are summarized from MODIS thermal anomaly and fire detection product (MOD14A1) and Burned Area Product(MCD64)',
                'MODIS MOD14A1 Data Source': 'https://modis.gsfc.nasa.gov/data/dataprod/mod14.php', 
                'Original MODIS MCD64 Data Source': 'https://modis.gsfc.nasa.gov/data/dataprod/mod45.php'
                }
    fds.to_netcdf(os.path.join(indir,'ComExDBM/04_OutputData','02_FireData','ComExDBM_'+str(y)+'_Fires'+'_V01.nc'))
    #
    dds = ydata[dcols].set_index(['date', 'lon', 'lat']).to_xarray()
    dds.attrs = {'title': 'gridMET Drought Index', 'version': 'V1.0',
                'description': 'Gridded every 5-day drought index,including pdsi,spi,spei related index and category',
                'SpatialCoverage': 'contiguous United States', 'SpatialResolution': '0.5 x 0.5', 
                'TemporalRange': str(y)+'-01-01 -> '+str(y)+'-12-31',  'TemporalResolution': 'daily', 
                'Format': 'NetCDF-4/HDF-5', 
                'comment': 'the drought index are regrided based on data orininally from the gridMET product',
                'Original gridMET Data Source': 'https://www.climatologylab.org/gridmet.html'
                }
    dds.to_netcdf(os.path.join(indir,'ComExDBM/04_OutputData','03_DroughtData','ComExDBM_'+str(y)+'_Drought'+'_V01.nc'))
    #
    ndata = pd.read_csv(os.path.join(indir,'NARR_vars_daily','NARR_regrid','NARR_vars_'+str(y)+'.csv'))
    ncols = ['date', 'lon', 'lat']+[c for c in ndata.columns if any(k in c for k in ['air_sfc','soilm', 'lhtfl', 'shtfl','apcp'])]
    #
    mdata = pd.read_csv(os.path.join(indir,'MERRA2_vars_daily','MERRA_regrid','MERRA_asm_'+str(y)+'.csv'))
    mcols = ['date', 'lon', 'lat']+[c for c in mdata.columns if any(k in c for k in ['RH','QV','U','V','BCOC'])]
    # merge NARR+MERRA2
    mndata = ndata[ncols].merge(mdata[mcols],on= ['date', 'lon', 'lat'],how='outer')
    mnds = mndata.set_index(['date', 'lon', 'lat']).to_xarray()
    mnds.attrs= {'title': 'Meteorological variables from NARR and MERRA2','version': 'V1.0', 
                'NARR description': 'daily summary(min/max/mean) of NARR 3-hourly Meteorological variables(i.e.,air_sfc->air temperature,soilm->soil misture,lhtfl->latent heat flux, shtfl->sensible heat flux, apcp->accumulated precipitation )',  
                'MERRA2 description': 'daily summary(min/max/mean) of MERRA2 3-hourly Meteorological variables(i.e.,air_sfc->air temperature,soilm->soil misture,lhtfl->latent heat flux, shtfl->sensible heat flux, apcp->accumulated precipitation )',             
                'SpatialCoverage': 'contiguous United States', 'SpatialResolution': '0.5 x 0.5', 
                'TemporalRange': str(y)+'-01-01 -> '+str(y)+'-12-31',  'TemporalResolution': 'daily', 
                'Format': 'NetCDF-4/HDF-5', 
                'comment': 'the original 3-hourly data are downloaded from ',
                'NARR references': 'https://www.esrl.noaa.gov/psd/data/gridded/data.narr.html', 'MERRA2 References': 'http://gmao.gsfc.nasa.gov',
                'Original NARR data source': 'http://www.emc.ncep.noaa.gov/mmb/rreanl/index.html',
                'Original MERRA2 data source': 'https://disc.gsfc.nasa.gov/datasets/M2I3NPASM_5.12.4/summary'}
    mnds.to_netcdf(os.path.join(indir,'ComExDBM/04_OutputData','04_Meteorological_Variables','ComExDBM_'+str(y)+'_MetVars'+'_V01.nc'))





# nds = ndata[ncols].set_index(['date', 'lon', 'lat']).to_xarray()
# nds.attrs= {'title': 'NARR Meteorological variables','version': 'V1.0', 
            # 'dataset_title': 'NCEP North American Regional Reanalysis (NARR)', 
            # 'description': 'daily summary(min/max/mean) of NARR 3-hourly Meteorological variables(i.e.,air_sfc->air temperature,soilm->soil misture,lhtfl->latent heat flux, shtfl->sensible heat flux, apcp->accumulated precipitation )',      
            # 'SpatialCoverage': 'contiguous United States', 'SpatialResolution': '0.5 x 0.5', 
            # 'TemporalRange': str(y)+'-01-01 -> '+str(y)+'-12-31',  'TemporalResolution': 'daily', 
            # 'Format': 'NetCDF-4/HDF-5', 
            # 'comment': 'the original 3-hourly data are created Mon Jul 18 13:06:52 MDT 2016 by NOAA/ESRL/PSD',
            # 'references': 'https://www.esrl.noaa.gov/psd/data/gridded/data.narr.html', 
            # 'Original data source': 'http://www.emc.ncep.noaa.gov/mmb/rreanl/index.html', 'References': ''}
# nds.to_netcdf(os.path.join(indir,'ComExDBM/04_OutputData','04_NARR','ComExDBM_'+str(y)+'_NARR_asm'+'_V01.nc'))
# #
# mdata = pd.read_csv(os.path.join(indir,'MERRA2_vars_daily','MERRA_regrid','MERRA_asm_'+str(y)+'.csv'))
# mcols = ['date', 'lon', 'lat']+[c for c in mdata.columns if any(k in c for k in ['RH','QV','U','V','BCOC'])]
# mds = mdata[mcols].set_index(['date', 'lon', 'lat']).to_xarray()
# mds.attrs= {'Title': 'MERRA2 Meteorological variables', 'version': 'V1.0', 
            # 'dataset_title': 'Modern-Era Retrospective analysis for Research and Applications version 2 (MERRA-2)', 
            # 'description': 'daily summary(min/max/mean) of MERRA2 3-hourly Meteorological variables(i.e.,RH->relative humidity, QV->specific humidity, U->Uwind,V->Vwind,BCOC->carbon aerosol)',      
            # 'comment': 'the original file generated: Sun May 13 22:19:23 2018 GMT', 
            # 'SpatialCoverage': 'contiguous United States', 'SpatialResolution': '0.5 x 0.5', 
            # 'TemporalRange': str(y)+'-01-01 -> '+str(y)+'-12-31',  'TemporalResolution': 'daily', 
            # 'Format': 'NetCDF-4/HDF-5', 
            # 'References': 'http://gmao.gsfc.nasa.gov',
            # 'Original data source': 'http://www.emc.ncep.noaa.gov/mmb/rreanl/index.html', 'References': ''} 
# mds.to_netcdf(os.path.join(indir,'ComExDBM/04_OutputData','05_MERRA2','ComExDBM_'+str(y)+'_MERRA2_asm'+'_V01.nc'))