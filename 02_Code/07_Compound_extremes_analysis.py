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
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.ticker as ticker
from matplotlib.gridspec import GridSpec
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
# import metpy.calc as mpcalc
# from metpy.units import units
import itertools

from copulas.multivariate import GaussianMultivariate,VineCopula
from sklearn.preprocessing import MinMaxScaler

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


##############################################################################
if __name__ == "__main__":
    ##PNW(WA,ID,OR)
# file system
indir = '/qfs/projects/comexdbm/LDRD_extremes/'
conus = gpd.read_file(os.path.join(indir,'CONUS4326','contig_us.shp'), crs="EPSG:4326")
years= list(range(2001, 2021))
months =[str(i).zfill(2) for i in list(range(1, 13))]
PNW_shape = conus[conus.NAME.isin(['Washington','Idaho','Oregon'])]
PW_shape = conus[conus.NAME.isin(['Washington','Idaho','Oregon','California'])]
SEUS_shape = conus[conus.NAME.isin(['Texas','Oklahoma','Arkansas','Louisiana','Mississippi','Alabama'])]
    
years= list(range(2001, 2021))
ydf = []
for y in years:
    #ydata = tdf[tdf.date.dt.year==y]
    #fdata['discovery_date']=pd.to_datedate(fdata['discovery_date']);
    #yfile = 'CPC_Global_Unified_Temperature_with_HW_labels_Fires_Drought_'+str(y)+'.csv'
    yfile = 'CPC_Global_Unified_Temperature_HFD_'+str(y)+'.csv'
    ydata = pd.read_csv(os.path.join(indir,'HW_Wildfire_Drought_data',yfile))#.iloc[:, 1:]
    ydata['date']=pd.to_datetime(ydata['date']);
    #
    hindex = ['HI02','HI04','HI05','HI06','HI09','HI10']
    findex= ['FHS_c8c9','FHS_c9']
    mindex =['maxFRP_max_c8c9','maxFRP_max_c9']
    #sdata= ydata[ydata[hindex].values==1]
    #sdata=ydata[(ydata.date.dt.month>=5)&((ydata.date.dt.month<=10))]
    sdata=ydata
    sdata =sdata.dropna(subset =['tmax', 'tmin'])
    sdata['year'] = sdata.date.dt.year
    sdata['month'] = sdata.date.dt.month
    ydf+=[sdata.copy()]
    # sdf1 =sdata.groupby(['lon', 'lat','year'])[hindex].sum().reset_index()
    # sdf2 =sdata.groupby(['lon', 'lat','year'])[findex].sum().reset_index()
    # sdf3 =sdata.groupby(['lon', 'lat','year'])[mindex].max().reset_index()
    # #sdf=sdf1.merge(sdf2,on=['lon', 'lat','year'])
    # sdf=reduce(lambda x, y: pd.merge(x, y, on = ['lat', 'lon','year']), [sdf1,sdf2,sdf3])
    #ydf+=[sdata.copy()]

##
ydf=pd.concat(ydf)
ydf['FHS_c9_cat'] =0
ydf['FHS_c9_cat'][ydf.FHS_c9>0]=1


hindex = ['HI02','HI04','HI05','HI06','HI09','HI10']
findex= ['FHS_c9','FHS_c8c9']
mindex =['maxFRP_max_c9']
dindex = ['pdsi','spi14d','spi30d','spi90d','spei14d','spei30d','spei90d']
dindex1 = ['pdsi','spi14d','spi30d','spi90d']
dindex2 = ['spei14d','spei30d','spei90d']
cindex = [h +'_FHS_c9' for h in hindex]

for d in dindex1+dindex2:
    ydf[d]=0
    ydf[d][ydf[d+'_category']<=1]=1

for h in hindex : 
    ydf[h+'_FHS_c9'] =0
    ydf[h+'_FHS_c9'][(ydf[h]==1)&(ydf.FHS_c9==0)] =1
    ydf[h+'_FHS_c9'][(ydf[h]==0)&(ydf.FHS_c9>0)] =2
    ydf[h+'_FHS_c9'][(ydf[h]==1)&(ydf.FHS_c9>0)] =3

for h in dindex1+dindex2:
    ydf[h+'_FHS_c9'] =0
    ydf[h+'_FHS_c9'][(ydf[h]==1)&(ydf.FHS_c9==0)] =1
    ydf[h+'_FHS_c9'][(ydf[h]==0)&(ydf.FHS_c9>0)] =2
    ydf[h+'_FHS_c9'][(ydf[h]==1)&(ydf.FHS_c9>0)] =3


###########################
#sdf = ydf[(ydf.time.dt.year==2018)&(ydf.time.dt.month==7)]
#23-29 July 2018
sdf = ydf[(ydf['time'] > '2006-07-12') & (ydf['time'] < '2006-07-19')]
cps = sdf.groupby(['lon','lat','time']).agg({i:'value_counts' for i in cindex}).reset_index()
#sdates = ['2018-07-23':'2018-07-30']
#for cp in cindex :
cps.time= cps.time.dt.strftime('%Y-%m-%d')
sdates = cps.time.unique()
cp= cindex[0]
for s in sdates:
    cdata = cps[cps.time==s][['lon','lat','level_3',cp]]
    data=cdata.dropna()
    #data =cdata[cdata.level_2==3]
    fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(6, 6))
    #ax = conus_new.plot(facecolor="none",edgecolor=conus['region'],linewidth=2)
    #ax.set_title('Correlation between '+h+' days'+' and '+ 'Fire Grids '+f,fontsize=8)
    #divider = make_axes_locatable(ax)
    #cax = divider.append_axes("right", size="3%", pad=0.01)
    x, y = data['lon'].values, data['lat'].values
    # x, y = tdata_ngp['lon'].values, tdata_ngp['lat'].values
    # im=ax.scatter(-123.75, 46.75, c='#EE3B3B', marker="*",s=25,alpha=0.6) #data[(data.lon==-122.75)&(data.lat==47.75)]['corr'],marker="s",alpha=0.8, cmap=color_ramp
    # im=ax.scatter(-120.75, 45.25, c='#EE3B3B', marker="*",s=25,alpha=0.6) #data[(data.lon==-117.75)&(data.lat==44.75)]['corr']
    # im=ax.scatter(-115.25, 46.25, c='#EE3B3B', marker="*",s=25,alpha=0.6) #data[(data.lon==-111.75)&(data.lat==42.25)]['corr']
    colors = {0:'white', 1:'skyblue', 2:'orange', 3:'red'}
    im=ax.scatter(x, y, c=data['level_3'].map(colors), marker="s",s=10,alpha=0.8) #,marker="s",alpha=0.8, cmap=color_ramp #cmap=plt.get_cmap("OrRd") #
    # add a legend
    handles = [Line2D([0], [0], marker='s', color='w', markerfacecolor=v, label=k, markersize=6) for k, v in colors.items()]
    ax.legend(title='Extremes', handles=handles,labels=['None', 'H', 'F','HF'],  loc='lower right')
    #im=ax.scatter(x, y, marker="*",s=0.01) #,marker="s"
    #im.set_clim(data['gl'].quantile(.1),data['gl'].quantile(.9))
    #im.set_clim(0,4)
    conus.plot(facecolor="none",edgecolor='r',linewidth=1,alpha=0.2,ax=ax)
    #PNW_shape.plot(facecolor="none",edgecolor='k',linewidth=1,ax=ax)
    #fig.colorbar(im, cax=cax, orientation='vertical',shrink=0.9)
    ax.set_xlabel('Lon');ax.set_ylabel('Lat')
    plt.title('Compound extremes'+'('+s+')')
    #plt.legend(handles=handles, title='Extremes')
    plt.savefig(os.path.join('/home/linx882/LDRD_extremes/','fire_plots','compound_'+cp+'_'+s+'.png'),
                transparent = False,facecolor = 'white',bbox_inches='tight',dpi=300)

######################################################################
from scipy.stats import boxcox
import numpy as np
#from copulas.bivariate.gumbel import Gumbel
from copulas.bivariate import Gumbel
from scipy.stats import kendalltau
from copulas.visualization import scatter_3d
from scipy.stats import rankdata
from copulas.bivariate import Frank, Clayton

cdata =ydf[(ydf['tmax']>34.6)&(ydf['FHS_c9']>=21)] #ydf[(ydf['HI02_FHS_c9']!=0)&(ydf['FHS_c9']>=59)] #~95% &(ydf['FHS_c9']>=38)&(ydf['tmax']>=59)
data = cdata[['tmax','FHS_c9']]
scaler = MinMaxScaler()
normalized_data = scaler.fit_transform(data)
normalized_data = pd.DataFrame(normalized_data, columns=data.columns)

# Compute the Kendall's tau
tau, _ = kendalltau(normalized_data['tmax'], normalized_data['FHS_c9'])

# Convert Kendall's tau to the Gumbel copula parameter
theta = 1 + 1 / (1 - tau)

# Fit the Gumbel copula
fc=  Frank()# Gumbel()
# gumbel.tau = tau
# gumbel.theta = theta
fc.fit(normalized_data[['tmax','FHS_c9']].values)

#u = gumbel.probability(normalized_data[['tmax','FHS_c9']].values)
u = fc.cdf(normalized_data[['tmax','FHS_c9']].values)
u_rank = np.apply_along_axis(rankdata, 1, u.reshape(-1, 1))
utd = np.corrcoef(u_rank[:, 0], u_rank[:, 1])[0, 1]
print("Upper tail dependence coefficient:", utd)

print("Lower tail dependence:", lower_tail_dependence)
print("Upper tail dependence:", upper_tail_dependence)


cdata =ydf[(ydf['tmax']>34.6)&(ydf['FHS_c9']>=21)] #ydf[(ydf['HI02_FHS_c9']!=0)&(ydf['FHS_c9']>=59)] #~95% &(ydf['FHS_c9']>=38)&(ydf['tmax']>=59)
data = cdata[['tmax','FHS_c9']]
scaler = MinMaxScaler()
normalized_data = scaler.fit_transform(data)
normalized_data = pd.DataFrame(normalized_data, columns=data.columns)


data['tt'], _ = boxcox(data.tmax)
data['ft'], _ = boxcox(data.FHS_c9)

copula = GaussianMultivariate()#VineCopula('regular')#
# copula.fit(normalized_data)
# joint_cdf = copula.probability_density(normalized_data)
copula.fit(data[['tt','ft']])
joint_cdf = copula.probability_density(data[['tt','ft']])

plt.figure(figsize=(8, 6))
#plt.scatter(data['FHS_c9'],data['tmax'], c=joint_cdf, cmap='viridis', marker='o', edgecolors='k', s=50)
plt.scatter(data['tt'],data['ft'], c=joint_cdf, cmap='viridis', marker='o', edgecolors='k', s=50)
plt.ylabel('Tmax')
plt.xlabel('FHS')
plt.colorbar(label='Joint PDF')
plt.title('Compound Extreme Events (Tmax, Fires) with Gaussian Copula')
plt.show()


fig = plt.figure(figsize=(8, 6))
gs = GridSpec(4, 4)
ax_scatter = fig.add_subplot(gs[1:4, 0:3])
ax_hist_y = fig.add_subplot(gs[0,0:3])
ax_hist_x = fig.add_subplot(gs[1:4, 3])
ax_scatter.scatter(data['FHS_c9'],data['tmax'], marker='o', edgecolors='r', s=5)
ax_scatter.set_xlabel('FHS (>95%)',color='k')
ax_scatter.set_ylabel('Tmax(>95%)',color='k')
ax_hist_y.hist(data['FHS_c9'], color='tab:red', alpha=0.4)
ax_hist_x.hist(data['tmax'], orientation = 'horizontal', color='tab:red', alpha=0.4)
plt.savefig(os.path.join(cdir,'scatter_Tmax_FHS_c9'+'.png'),
            transparent = False,facecolor = 'white',bbox_inches='tight',dpi=300)
plt.show()

#data = cdata[['tmax','FHS_c9']]
data = cdata[['tmax','FHS_c9']]
#data['FHS_c9'] =np.log10(data['FHS_c9']+1);data['tmax'] =np.log10(data['tmax']+1)
fig = plt.figure(figsize=(8, 6))
gs = GridSpec(4, 4)
ax_scatter = fig.add_subplot(gs[1:4, 0:3])
ax_hist_y = fig.add_subplot(gs[0,0:3])
ax_hist_x = fig.add_subplot(gs[1:4, 3])
ax_scatter.scatter(data['tmax'],data['FHS_c9'], marker='o', edgecolors='r', s=5)
ax_scatter.set_ylabel('FHS(>95%) ',color='k')
ax_scatter.set_xlabel('Tmax(>95%)',color='k')
ax_hist_x.hist(data['FHS_c9'], color='tab:red', alpha=0.4)
ax_hist_y.hist(data['tmax'], orientation = 'horizontal', color='tab:red', alpha=0.4)
plt.savefig(os.path.join(cdir,'scatter_Tmax_FHS_c9_95'+'.png'),
            transparent = False,facecolor = 'white',bbox_inches='tight',dpi=300)
plt.show()

data = ydf[['tmax','FHS_c9']]
#data['FHS_c9'] =np.log10(data['FHS_c9']+1);data['tmax'] =np.log10(data['tmax']+1)
fig = plt.figure(figsize=(8, 6))
gs = GridSpec(4, 4)
ax_scatter = fig.add_subplot(gs[1:4, 0:3])
ax_hist_y = fig.add_subplot(gs[0,0:3])
ax_hist_x = fig.add_subplot(gs[1:4, 3])
ax_scatter.scatter(data['tmax'],data['FHS_c9'], marker='o', edgecolors='r', s=5)
ax_scatter.set_ylabel('FHS ',color='k')
ax_scatter.set_xlabel('Tmax',color='k')
ax_hist_x.hist(data['FHS_c9'], color='tab:red', alpha=0.4)
ax_hist_y.hist(data['tmax'], orientation = 'horizontal', color='tab:red', alpha=0.4)
plt.savefig(os.path.join(cdir,'scatter_Tmax_FHS_c9_95_all'+'.png'),
            transparent = False,facecolor = 'white',bbox_inches='tight',dpi=300)




sns.scatterplot(x='FHS_c9', y='tmax', data=ydf[(ydf['tmax']>34.6)&(ydf['FHS_c9']>=21)], alpha=0.5)
cdata = ydf[(ydf['HI02_FHS_c9']!=0)] 
sns.scatterplot(x='FHS_c9', y='tmax', data=cdata, alpha=0.5)

##############################
cdir= os.path.join(indir,'Cprob_plots')

for h in hindex : 
    fdata9 =ydf[ydf.FHS_c9>0]
    #fdata9.groupby(['column_1', 'column_2'])['target_column'].value_counts(normalize=True)
    cp1 = fdata9.groupby(['lon', 'lat'])[h].value_counts(normalize=True)
    cp1= cp1.unstack(fill_value=0).stack().reset_index()
    cp1=cp1[cp1[h]==1]
    data = cp1[cp1[0]>0]
    fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(5, 5))
    #ax = conus_new.plot(facecolor="none",edgecolor=conus['region'],linewidth=2)
    #ax.set_title('Correlation between '+h+' days'+' and '+ 'Fire Grids '+f,fontsize=8)
    # divider = make_axes_locatable(ax)
    # cax = divider.append_axes("right", size="3%", pad=0.01)
    x, y = data['lon'].values, data['lat'].values
    # x, y = tdata_ngp['lon'].values, tdata_ngp['lat'].values
    # im=ax.scatter(-123.75, 46.75, c='#EE3B3B', marker="*",s=25,alpha=0.6) #data[(data.lon==-122.75)&(data.lat==47.75)]['corr'],marker="s",alpha=0.8, cmap=color_ramp
    # im=ax.scatter(-120.75, 45.25, c='#EE3B3B', marker="*",s=25,alpha=0.6) #data[(data.lon==-117.75)&(data.lat==44.75)]['corr']
    # im=ax.scatter(-115.25, 46.25, c='#EE3B3B', marker="*",s=25,alpha=0.6) #data[(data.lon==-111.75)&(data.lat==42.25)]['corr']
    im=ax.scatter(x, y, c=data[0], marker="s",s=3,cmap=plt.get_cmap("OrRd"),alpha=0.8) #,marker="s",alpha=0.8, cmap=color_ramp
    #im=ax.scatter(x, y, marker="*",s=0.01) #,marker="s"
    #im.set_clim(data['gl'].quantile(.1),data['gl'].quantile(.9))
    im.set_clim(0,0.5)
    conus.plot(facecolor="none",edgecolor='k',linewidth=1,alpha=0.2,ax=ax)
    #fig.colorbar(im, cax=cax, orientation='vertical',shrink=0.9)
    axins = inset_axes(
                        ax,
                        width="5%",  # width: 5% of parent_bbox width
                        height="100%",  # height: 50%
                        loc="lower left",
                        bbox_to_anchor=(1.02, 0., 1, 1),
                        bbox_transform=ax.transAxes,
                        borderpad=0,
                        )
    cbar = fig.colorbar(im, cax=axins, orientation='vertical',shrink=0.9)
    cbar.ax.locator_params(nbins=5)
    cbar.formatter.set_scientific(True)
    cbar.formatter.set_powerlimits((0, 0))
    ax.set_xlabel('Lon')
    ax.set_ylabel('Lat')
    #PNW_shape.plot(facecolor="none",edgecolor='k',linewidth=1,ax=ax)
    plt.savefig(os.path.join(cdir,'cprob_'+h+'_FHS_c9_5'+'.png'),
                transparent = False,facecolor = 'white',bbox_inches='tight',dpi=300)




##############################
hdf = ydf[ydf.HI02==1]
cp1 = hdf.groupby(['lon', 'lat'])['FHS_c9_cat'].value_counts(normalize=True)
cp1= cp1.unstack(fill_value=0).stack().reset_index()
cp1=cp1[cp1.FHS_c9_cat==1]



for h in hindex : 
    fdata9 =ydf[ydf[h]==1]
    #fdata9.groupby(['column_1', 'column_2'])['target_column'].value_counts(normalize=True)
    cp1 = fdata9.groupby(['lon', 'lat'])['FHS_c9_cat'].value_counts(normalize=True)
    cp1= cp1.unstack(fill_value=0).stack().reset_index()
    cp1=cp1[cp1['FHS_c9_cat']==1]
    data = cp1[cp1[0]>0]
    fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(5, 5))
    #ax = conus_new.plot(facecolor="none",edgecolor=conus['region'],linewidth=2)
    #ax.set_title('Correlation between '+h+' days'+' and '+ 'Fire Grids '+f,fontsize=8)
    # divider = make_axes_locatable(ax)
    # cax = divider.append_axes("right", size="3%", pad=0.01)
    x, y = data['lon'].values, data['lat'].values
    # x, y = tdata_ngp['lon'].values, tdata_ngp['lat'].values
    # im=ax.scatter(-123.75, 46.75, c='#EE3B3B', marker="*",s=25,alpha=0.6) #data[(data.lon==-122.75)&(data.lat==47.75)]['corr'],marker="s",alpha=0.8, cmap=color_ramp
    # im=ax.scatter(-120.75, 45.25, c='#EE3B3B', marker="*",s=25,alpha=0.6) #data[(data.lon==-117.75)&(data.lat==44.75)]['corr']
    # im=ax.scatter(-115.25, 46.25, c='#EE3B3B', marker="*",s=25,alpha=0.6) #data[(data.lon==-111.75)&(data.lat==42.25)]['corr']
    im=ax.scatter(x, y, c=data[0], marker="s",s=3,cmap=plt.get_cmap("OrRd"),alpha=0.8) #,marker="s",alpha=0.8, cmap=color_ramp
    #im=ax.scatter(x, y, marker="*",s=0.01) #,marker="s"
    #im.set_clim(data['gl'].quantile(.1),data['gl'].quantile(.9))
    #im.set_clim(0,data[0].max()*0.75)
    #im.set_clim(0,0.2)
    if h=='HI04':
        im.set_clim(0,0.2)
    else:
        im.set_clim(0,0.1)
    conus.plot(facecolor="none",edgecolor='k',linewidth=1,alpha=0.2,ax=ax)
    #PNW_shape.plot(facecolor="none",edgecolor='k',linewidth=1,ax=ax)
    #fig.colorbar(im, cax=cax, orientation='vertical',shrink=0.9)
    axins = inset_axes(
                        ax,
                        width="5%",  # width: 5% of parent_bbox width
                        height="100%",  # height: 50%
                        loc="lower left",
                        bbox_to_anchor=(1.02, 0., 1, 1),
                        bbox_transform=ax.transAxes,
                        borderpad=0,
                        )
    cbar = fig.colorbar(im, cax=axins, orientation='vertical',shrink=0.9)
    cbar.ax.locator_params(nbins=5)
    cbar.formatter.set_scientific(True)
    cbar.formatter.set_powerlimits((0, 0))
    ax.set_xlabel('Lon')
    ax.set_ylabel('Lat')
    plt.savefig(os.path.join(cdir,'cprob_FHS_c9_'+h+'.png'),
                transparent = False,facecolor = 'white',bbox_inches='tight',dpi=300)


for h in hindex : 
    #fdata9 =ydf[ydf[h]==1]
    fdata9 =ydf #[ydf[h]==1]
    #fdata9.groupby(['column_1', 'column_2'])['target_column'].value_counts(normalize=True)
    # cp1 = fdata9.groupby(['lon', 'lat'])['FHS_c9_cat'].value_counts(normalize=True)
    # cp1= cp1.unstack(fill_value=0).stack().reset_index()
    # cp1=cp1[cp1['FHS_c9_cat']==1]
    cp1 = fdata9.groupby(['lon', 'lat'])[h].value_counts(normalize=True)
    cp1= cp1.unstack(fill_value=0).stack().reset_index()
    cp1=cp1[cp1[h]==1]
    data = cp1[cp1[0]>0]
    fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(5, 5))
    #ax = conus_new.plot(facecolor="none",edgecolor=conus['region'],linewidth=2)
    #ax.set_title('Correlation between '+h+' days'+' and '+ 'Fire Grids '+f,fontsize=8)
    # divider = make_axes_locatable(ax)
    # cax = divider.append_axes("right", size="3%", pad=0.01)
    x, y = data['lon'].values, data['lat'].values
    # x, y = tdata_ngp['lon'].values, tdata_ngp['lat'].values
    # im=ax.scatter(-123.75, 46.75, c='#EE3B3B', marker="*",s=25,alpha=0.6) #data[(data.lon==-122.75)&(data.lat==47.75)]['corr'],marker="s",alpha=0.8, cmap=color_ramp
    # im=ax.scatter(-120.75, 45.25, c='#EE3B3B', marker="*",s=25,alpha=0.6) #data[(data.lon==-117.75)&(data.lat==44.75)]['corr']
    # im=ax.scatter(-115.25, 46.25, c='#EE3B3B', marker="*",s=25,alpha=0.6) #data[(data.lon==-111.75)&(data.lat==42.25)]['corr']
    im=ax.scatter(x, y, c=data[0], marker="s",s=3,cmap=plt.get_cmap("OrRd"),alpha=0.6) #,marker="s",alpha=0.8, cmap=color_ramp
    #im=ax.scatter(x, y, marker="*",s=0.01) #,marker="s"
    #im.set_clim(data['gl'].quantile(.1),data['gl'].quantile(.9))
    #im.set_clim(0,data[0].max())
    #im.set_clim(0,0.1)
    if h=='HI04':
        im.set_clim(0,0.01)
    elif h=='HI10':
        im.set_clim(0,0.1)
    else:
        im.set_clim(0,0.06)
    conus.plot(facecolor="none",edgecolor='k',linewidth=1,alpha=0.2,ax=ax)
    #PNW_shape.plot(facecolor="none",edgecolor='k',linewidth=1,ax=ax)
    #fig.colorbar(im, cax=cax, orientation='vertical',shrink=0.9)
    axins = inset_axes(
                        ax,
                        width="5%",  # width: 5% of parent_bbox width
                        height="100%",  # height: 50%
                        loc="lower left",
                        bbox_to_anchor=(1.02, 0., 1, 1),
                        bbox_transform=ax.transAxes,
                        borderpad=0,
                        )
    cbar = fig.colorbar(im, cax=axins, orientation='vertical',shrink=0.9)
    cbar.ax.locator_params(nbins=5)
    cbar.formatter.set_scientific(True)
    cbar.formatter.set_powerlimits((0, 0))
    ax.set_xlabel('Lon')
    ax.set_ylabel('Lat')
    plt.savefig(os.path.join(cdir,'Mprob_'+h+'.png'),
                transparent = False,facecolor = 'white',bbox_inches='tight',dpi=300)



fdata9 =ydf #[ydf[h]==1]
#fdata9.groupby(['column_1', 'column_2'])['target_column'].value_counts(normalize=True)
cp1 = fdata9.groupby(['lon', 'lat'])['FHS_c9_cat'].value_counts(normalize=True)
cp1= cp1.unstack(fill_value=0).stack().reset_index()
cp1=cp1[cp1['FHS_c9_cat']==1]
# cp1 = fdata9.groupby(['lon', 'lat'])[h].value_counts(normalize=True)
# cp1= cp1.unstack(fill_value=0).stack().reset_index()
# cp1=cp1[cp1[h]==1]
data = cp1[cp1[0]>0]
fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(5, 5))
#ax = conus_new.plot(facecolor="none",edgecolor=conus['region'],linewidth=2)
#ax.set_title('Correlation between '+h+' days'+' and '+ 'Fire Grids '+f,fontsize=8)
# divider = make_axes_locatable(ax)
# cax = divider.append_axes("right", size="3%", pad=0.01)
x, y = data['lon'].values, data['lat'].values
# x, y = tdata_ngp['lon'].values, tdata_ngp['lat'].values
# im=ax.scatter(-123.75, 46.75, c='#EE3B3B', marker="*",s=25,alpha=0.6) #data[(data.lon==-122.75)&(data.lat==47.75)]['corr'],marker="s",alpha=0.8, cmap=color_ramp
# im=ax.scatter(-120.75, 45.25, c='#EE3B3B', marker="*",s=25,alpha=0.6) #data[(data.lon==-117.75)&(data.lat==44.75)]['corr']
# im=ax.scatter(-115.25, 46.25, c='#EE3B3B', marker="*",s=25,alpha=0.6) #data[(data.lon==-111.75)&(data.lat==42.25)]['corr']
im=ax.scatter(x, y, c=data[0], marker="s",s=3,cmap=plt.get_cmap("OrRd"),alpha=0.6) #,marker="s",alpha=0.8, cmap=color_ramp
#im=ax.scatter(x, y, marker="*",s=0.01) #,marker="s"
#im.set_clim(data['gl'].quantile(.1),data['gl'].quantile(.9))
#im.set_clim(0,data[0].max())
im.set_clim(0,0.05)
conus.plot(facecolor="none",edgecolor='k',linewidth=1,alpha=0.2,ax=ax)
#cbar = fig.colorbar(im, cax=cax, orientation='vertical',shrink=0.9,pad = 0)
#PNW_shape.plot(facecolor="none",edgecolor='k',linewidth=1,ax=ax)
axins = inset_axes(
    ax,
    width="5%",  # width: 5% of parent_bbox width
    height="100%",  # height: 50%
    loc="lower left",
    bbox_to_anchor=(1.02, 0., 1, 1),
    bbox_transform=ax.transAxes,
    borderpad=0,
)
cbar = fig.colorbar(im, cax=axins, orientation='vertical',shrink=0.9)
cbar.ax.locator_params(nbins=5)
cbar.formatter.set_scientific(True)
cbar.formatter.set_powerlimits((0, 0))
ax.set_xlabel('Lon')
ax.set_ylabel('Lat')
plt.savefig(os.path.join(cdir,'Mprob_'+'FHS_c9_5'+'.png'),
            transparent = False,facecolor = 'white',bbox_inches='tight',dpi=300)



for h in hindex : 
    #fdata9 =ydf[ydf[h]==1]
    fdata9 =ydf #[ydf[h]==1]
    fdata9['fh'] =0
    fdata9['fh'] [(fdata9[h]==1)&(fdata9['FHS_c9_cat']==1)] =1
    #fdata9.groupby(['column_1', 'column_2'])['target_column'].value_counts(normalize=True)
    # cp1 = fdata9.groupby(['lon', 'lat'])['FHS_c9_cat'].value_counts(normalize=True)
    # cp1= cp1.unstack(fill_value=0).stack().reset_index()
    # cp1=cp1[cp1['FHS_c9_cat']==1]
    cp1 = fdata9.groupby(['lon', 'lat'])['fh'].value_counts(normalize=True)
    cp1= cp1.unstack(fill_value=0).stack().reset_index()
    cp1=cp1[cp1['fh']==1]
    data = cp1[cp1[0]>0]
    fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(5, 5))
    #ax = conus_new.plot(facecolor="none",edgecolor=conus['region'],linewidth=2)
    #ax.set_title('Correlation between '+h+' days'+' and '+ 'Fire Grids '+f,fontsize=8)
    # divider = make_axes_locatable(ax)
    # cax = divider.append_axes("right", size="3%", pad=0.01)
    x, y = data['lon'].values, data['lat'].values
    # x, y = tdata_ngp['lon'].values, tdata_ngp['lat'].values
    # im=ax.scatter(-123.75, 46.75, c='#EE3B3B', marker="*",s=25,alpha=0.6) #data[(data.lon==-122.75)&(data.lat==47.75)]['corr'],marker="s",alpha=0.8, cmap=color_ramp
    # im=ax.scatter(-120.75, 45.25, c='#EE3B3B', marker="*",s=25,alpha=0.6) #data[(data.lon==-117.75)&(data.lat==44.75)]['corr']
    # im=ax.scatter(-115.25, 46.25, c='#EE3B3B', marker="*",s=25,alpha=0.6) #data[(data.lon==-111.75)&(data.lat==42.25)]['corr']
    im=ax.scatter(x, y, c=data[0], marker="s",s=3,cmap=plt.get_cmap("OrRd"),alpha=1.0) #,marker="s",alpha=0.8, cmap=color_ramp
    #im=ax.scatter(x, y, marker="*",s=0.01) #,marker="s"
    #im.set_clim(data['gl'].quantile(.1),data['gl'].quantile(.9))
    #im.set_clim(0,0.1)
    if h=='HI04':
        im.set_clim(0,0.001)
    else:
        im.set_clim(0,0.005)
    conus.plot(facecolor="none",edgecolor='k',linewidth=1,alpha=0.2,ax=ax)
    #cbar = fig.colorbar(im, cax=cax, orientation='vertical',shrink=0.9,pad = 0)
    #PNW_shape.plot(facecolor="none",edgecolor='k',linewidth=1,ax=ax)
    axins = inset_axes(
        ax,
        width="5%",  # width: 5% of parent_bbox width
        height="100%",  # height: 50%
        loc="lower left",
        bbox_to_anchor=(1.02, 0., 1, 1),
        bbox_transform=ax.transAxes,
        borderpad=0,
    )
    cbar = fig.colorbar(im, cax=axins, orientation='vertical',shrink=0.9)
    cbar.ax.locator_params(nbins=5)
    cbar.formatter.set_scientific(True)
    cbar.formatter.set_powerlimits((0, 0))
    ax.set_xlabel('Lon')
    ax.set_ylabel('Lat')
    plt.savefig(os.path.join(cdir,'prob_'+h+'&FHS_c9_01'+'.png'),
                transparent = False,facecolor = 'white',bbox_inches='tight',dpi=300)




fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(5, 5))
x, y = data['lon'].values, data['lat'].values
im=ax.scatter(x, y,s=0.5)
conus.plot(facecolor="none",edgecolor='k',linewidth=1,alpha=0.2,ax=ax)
plt.show()

####################################################
fdf =ydf.copy()
#fdf =ydf[(ydf.state.isin(PNW))] 
hindex = ['HI02','HI04','HI05','HI06','HI09','HI10']
findex= ['FHS_c9','FHS_c8c9','BA_km2']
mindex =['maxFRP_max_c9']
dindex = ['pdsi','spi14d','spi30d','spi90d','spei14d','spei30d','spei90d']
dindex1 = ['pdsi','spi14d','spi30d','spi90d']
dindex2 = ['spei14d','spei30d','spei90d']
cindex = [h +'_FHS_c9' for h in hindex]


df1 = fdf.groupby(['year','lon','lat'])[hindex+dindex].sum().reset_index().groupby(['year'])[hindex+dindex].mean().reset_index()
#df2 = ydf.groupby(['year'])[mindex].max().reset_index()
df2 = fdf.groupby(['year','lon','lat'])[findex].sum().reset_index().groupby(['year'])[findex].mean().reset_index()
df3 = fdf.groupby(['year'])['tmax'].max().reset_index()
df4 = fdf.groupby(['year','lon','lat']).agg({i:'value_counts' for i in cindex}).reset_index()
df4 =df4[df4.level_3==3].groupby(['year'])[cindex].sum().reset_index()
df=reduce(lambda x, y: pd.merge(x, y, on = ['year']), [df1,df2,df3,df4])
df.set_index('year', inplace=True)

#fig, axes = plt.subplots(nrows=3, ncols=1,figsize=(6, 10))
fig, axes = plt.subplots(nrows=3, ncols=1,figsize=(6, 10))
ticklabels = ['']*len(df.index)
# Every 4th ticklable shows the month and day
#ticklabels[::6] = [item.strftime('%b %d') for item in df.index[::6]]
#df['tmax'].plot(ax=axes[0]); 
#axes[0].xaxis.set_major_formatter(ticker.FixedFormatter(ticklabels))
df[hindex].plot(ax=axes[0]);
axes[0].xaxis.set_ticks(np.arange(2001, 2021, 2))
axes[0].legend(ncol=2);
axes[0].set_xlabel('Year');axes[0].set_ylabel('Gridded annual average HW days')
#axes[1].xaxis.set_major_formatter(ticker.FixedFormatter(ticklabels))
df[findex[0:2]].plot(ax=axes[1]);#df[mindex].plot(ax=axes[1],secondary_y=True)
axes[1].xaxis.set_ticks(np.arange(2001, 2021, 2))
df['BA_km2'].plot(ax=axes[1],label='Burned Area(km2)', style='r--', secondary_y=True)
lines = axes[1].get_lines() + axes[1].right_ax.get_lines()
axes[1].legend(lines, [l.get_label() for l in lines],loc ='upper left', ncol=3) #,loc ='upper left'
axes[1].set_xlabel('Year');axes[1].set_ylabel('Gridded annual average FHS (BA) ')
#lines = axes[1].get_lines() + axes[1].right_ax.get_lines()
#axes[1].legend(lines, [l.get_label() for l in lines],loc ='upper left')
#axes[2].xaxis.set_major_formatter(ticker.FixedFormatter(ticklabels))
df[['spi14d','spi30d','spi90d']+dindex2].plot(ax=axes[2]); 
df['pdsi'].plot(ax=axes[2],label='pdsi', style='r--', secondary_y=True)
lines = axes[2].get_lines() + axes[2].right_ax.get_lines()
axes[2].legend(lines, [l.get_label() for l in lines],loc ='upper left', ncol=2) #,loc ='upper left'
axes[2].xaxis.set_ticks(np.arange(2001, 2021, 2))
axes[2].set_xlabel('Year');axes[2].set_ylabel('Gridded annual average Drought days')
# df[cindex].plot(ax=axes[2]);
# axes[2].xaxis.set_ticks(np.arange(2001, 2021, 2))
# axes[2].legend(ncol=3);
# axes[2].set_xlabel('Year');axes[2].set_ylabel('Annual sum of compound events')
# axes[2].set_ylim(0, 900)
#df[mindex].plot(ax=axes[3,0])
#df[dindex2].plot(ax=axes[3,1])
#axes[3].xaxis.set_major_formatter(ticker.FixedFormatter(ticklabels))
plt.tight_layout()
plt.savefig(os.path.join('/qfs/projects/comexdbm/LDRD_extremes/','fire_plots','Drought_HW_FireBA_ts_annual_average_conus_compound.png'),
                transparent = False,facecolor = 'white',bbox_inches='tight',dpi=300)
plt.show()


sdf= ydf[(ydf.FHS_c8c9>0)|(ydf.BA_km2>0)]
#fmax = sdf.groupby(['lat', 'lon'])[['FHS_c9','FHS_c8c9','BA_km2']].max().reset_index()
fmax = sdf.groupby(['lat', 'lon','year','month'])[['FHS_c9','FHS_c8c9','BA_km2']].sum().reset_index()
plt.plot(fmax.FHS_c8c9,fmax.BA_km2,'*')
#fmax = sdf.groupby(['lat', 'lon','date'])[['FHS_c9','FHS_c8c9','BA_km2']].sum().reset_index()

sdf= fmax
corrs = sdf.groupby(['lat', 'lon'])[['FHS_c9','BA_km2']].corr()
corrs =corrs.iloc[0::2][['BA_km2']]
corrs.reset_index(inplace=True);
corrs.columns = corrs.columns.str.replace('BA_km2', "corr") 

#
data= corrs
fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(5, 5))
x, y = data['lon'].values, data['lat'].values
im=ax.scatter(x, y, c=data['corr'], marker="s",s=3,cmap=plt.get_cmap("OrRd"),alpha=1.0) #,marker="s",alpha=0.8, cmap=color_ramp
im.set_clim(0,1)
conus.plot(facecolor="none",edgecolor='k',linewidth=1,alpha=0.2,ax=ax)
axins = inset_axes(
    ax,
    width="5%",  # width: 5% of parent_bbox width
    height="100%",  # height: 50%
    loc="lower left",
    bbox_to_anchor=(1.02, 0., 1, 1),
    bbox_transform=ax.transAxes,
    borderpad=0,
)
cbar = fig.colorbar(im, cax=axins, orientation='vertical',shrink=0.9)
cbar.ax.locator_params(nbins=5)
cbar.formatter.set_scientific(True)
cbar.formatter.set_powerlimits((0, 0))
ax.set_xlabel('Lon')
ax.set_ylabel('Lat')

plt.savefig(os.path.join('/qfs/projects/comexdbm/LDRD_extremes/','fire_plots','corr_BA&FHS_c9_monthly_sum'+'.png'),
            transparent = False,facecolor = 'white',bbox_inches='tight',dpi=300)
            
            
dsum = ydf.groupby(['date'])[['FHS_c9','FHS_c8c9','BA_km2']].sum().reset_index()
dsum.set_index('date', inplace=True)


fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(6, 4))
sns.regplot(x=dsum.FHS_c8c9, y=dsum.BA_km2,ci = None,marker='*',scatter_kws={'s':2},label='FHS_c8c9', color='blue',ax=ax)
sns.regplot(x=dsum.FHS_c9, y=dsum.BA_km2,ci = None,marker='+',scatter_kws={'s':2},label='FHS_c9', color='red',ax=ax)
plt.legend()
plt.xlabel('Daily total FHS over CONUS')
plt.ylabel('Daily total burned area(km2) over CONUS')
# plt.plot(dsum.FHS_c8c9,dsum.BA_km2,'*',label='FHS_c8c9')
# plt.plot(dsum.FHS_c9,dsum.BA_km2,'+',label='FHS_c9')
plt.savefig(os.path.join('/qfs/projects/comexdbm/LDRD_extremes/','fire_plots','scatter_BA&FHS_daily_sum'+'.png'),
            transparent = False,facecolor = 'white',bbox_inches='tight',dpi=300)
plt.show()

dsum[['FHS_c8c9','BA_km2']].corr()
dsum[['FHS_c8c9','BA_km2']][dsum.FHS_c8c9<3000].corr()


dsum[['FHS_c9','BA_km2']].corr()
dsum[['FHS_c9','BA_km2']][dsum.FHS_c9<2000].corr()


msum = ydf.groupby(['year','month'])[['FHS_c9','FHS_c8c9','BA_km2']].sum().reset_index()
#msum.set_index('date', inplace=True)

# fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(6, 4))
# #axes[1].xaxis.set_major_formatter(ticker.FixedFormatter(ticklabels))
# msum['FHS_c8c9'].plot(ax=ax);#df[mindex].plot(ax=axes[1],secondary_y=True)
# msum['BA_km2'].plot(ax=ax,label='Burned Area(km2)', style='r--', secondary_y=True)
# lines = ax.get_lines() + ax.right_ax.get_lines()
# ax.legend(lines, [l.get_label() for l in lines],loc ='upper left', ncol=2) #,loc ='upper left'
# ax.set_xlabel('Date');ax.set_ylabel('Daily total FHS (BA) ')
# plt.show()
fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(6, 4))
sns.regplot(x=msum.FHS_c8c9, y=msum.BA_km2,ci = None,marker='*',scatter_kws={'s':2},label='FHS_c8c9', color='blue',ax=ax)
sns.regplot(x=msum.FHS_c9, y=msum.BA_km2,ci = None,marker='+',scatter_kws={'s':2},label='FHS_c9', color='red',ax=ax)
plt.legend()
plt.xlabel('Monthly total FHS over CONUS')
plt.ylabel('Monthly total burned area(km2) over CONUS')
# plt.plot(dsum.FHS_c8c9,dsum.BA_km2,'*',label='FHS_c8c9')
# plt.plot(dsum.FHS_c9,dsum.BA_km2,'+',label='FHS_c9')
plt.savefig(os.path.join('/qfs/projects/comexdbm/LDRD_extremes/','fire_plots','scatter_BA&FHS_monthly_sum'+'.png'),
            transparent = False,facecolor = 'white',bbox_inches='tight',dpi=300)
plt.show()



msum[['FHS_c8c9','BA_km2']].corr()
msum[['FHS_c9','BA_km2']].corr()



fod = gpd.read_file(os.path.join(indir,'FPA_FOD_20210617','FPA_FOD_20210617.gpkg'))

fod.columns = [x.lower() for x in fod.columns]
# FIRE_SIZE = The estimate of acres within the final perimeter of the fire.
fod.fire_size = fod.fire_size/247.11
fod= fod[['fire_year','discovery_date', 'discovery_time',
       'fire_size', 'fire_size_class', 'latitude', 'longitude','state']]
fod.columns = ['year','date', 'time','fire_size', 'fire_size_class', 'latitude', 'longitude','state']
fod['date'] = pd.to_datetime(fod['date']).dt.tz_localize(None)
fod['month'] =fod.date.dt.month


fycount =fod[fod.year.isin(list(range(2001, 2019)))].groupby(['year'])['fire_size'].sum().reset_index()## rmse
fycount['ydays'] = fycount.fire_size/fycount.fire_size.sum()*100

fycount1 =ydf[ydf.year.isin(list(range(2001, 2019)))].groupby(['year'])['BA_km2'].sum().reset_index()## rmse
fycount1['ydays'] = fycount1.BA_km2/fycount1.BA_km2.sum()*100


fycount2 =ydf[ydf.year.isin(list(range(2001, 2019)))].groupby(['year'])['FHS_c9'].sum().reset_index()## rmse
fycount2.FHS_c9 = fycount2.FHS_c9*0.86
fycount2['ydays'] = fycount2.FHS_c9/fycount2.FHS_c9.sum()*100

fycount3 =ydf[ydf.year.isin(list(range(2001, 2019)))].groupby(['year'])['FHS_c8c9'].sum().reset_index()## rmse
fycount3.FHS_c8c9 = fycount3.FHS_c8c9*0.86
fycount3['ydays'] = fycount3.FHS_c8c9/fycount3.FHS_c8c9.sum()*100

# mfod['month'] =mfod.date.dt.month
# mfod['year'] =mfod.date.dt.year
# fycount4 =mfod[mfod.year.isin(list(range(2001, 2019)))].groupby(['year'])['daily_area_km2'].sum().reset_index()## rmse

fycount = fycount.merge(fycount1, on='year', how='left')
fycount = fycount.merge(fycount2, on='year', how='left')
fycount = fycount.merge(fycount3, on='year', how='left')

fmcount= fod[fod.year.isin(list(range(2001, 2019)))].groupby(['month'])['fire_size'].sum().reset_index()## rmse
fmcount['mdays'] = fmcount.fire_size/fmcount.fire_size.sum()*100

fmcount1 =ydf[ydf.year.isin(list(range(2001, 2019)))].groupby(['month'])['BA_km2'].sum().reset_index()## rmse
fmcount1['mdays'] = fmcount1.BA_km2/fmcount1.BA_km2.sum()*100

fmcount2 =ydf[ydf.year.isin(list(range(2001, 2019)))].groupby(['month'])['FHS_c9'].sum().reset_index()## rmse
fmcount2.FHS_c9 = fmcount2.FHS_c9*0.86
fmcount2['mdays'] = fmcount2.FHS_c9/fmcount2.FHS_c9.sum()*100

fmcount3 =ydf[ydf.year.isin(list(range(2001, 2019)))].groupby(['month'])['FHS_c8c9'].sum().reset_index()## rmse
fmcount3.FHS_c8c9 = fmcount3.FHS_c8c9*0.86
fmcount3['mdays'] = fmcount3.FHS_c8c9/fmcount3.FHS_c8c9.sum()*100


fmcount = fmcount.merge(fmcount1, on='month', how='left')
fmcount = fmcount.merge(fmcount2, on='month', how='left')
fmcount = fmcount.merge(fmcount3, on='month', how='left')
# # fmcount
# # fycount
#fmcount = fmcount[fmcount.month.isin([3,4,5,6,7,8,9])] 
# Plot fire and hail matches (add hail reports)
fig, axes = plt.subplots(figsize=(7,7), dpi=300, sharey=True)
ax1 = plt.subplot(2, 1, 1)
#ax.bar(gycount.year ,gycount.lon)
ax1.plot(fycount.year, fycount.fire_size,'-o',label='FOD')
#ax1.plot(fycount.year, fycount.FHS_c9,'-*',label='MOD14_c9')
#ax1.plot(fycount.year, fycount.FHS_c8c9,'-*',label='MOD14_c8c9')
ax1.plot(fycount.year, fycount.BA_km2,'-+',label='MCD64')
ax1.legend(loc="upper left",ncol=3)
ax1.set_xticks([2001,  2004, 2007, 2010,  2013,  2016])
ax1.set_xticklabels(labels=[2001,  2004, 2007, 2010,  2013,  2016])
ax1.set_ylim(10000, 45000);
ax1.set_xlabel('Year', fontsize=12)
ax1.set_ylabel('Total burned area (km^2)', fontsize=12)
ax2 = plt.subplot(2, 1, 2)
ax2.plot(fmcount.month, fmcount.fire_size,'-o',label='FOD')
#ax2.plot(fmcount.month, fmcount.FHS_c9,'-*',label='MOD14_c9')
#ax2.plot(fmcount.month, fmcount.FHS_c8c9,'-*',label='MOD14_c8c9')
ax2.plot(fmcount.month, fmcount.BA_km2,'-+',label='MCD64')
ax2.legend(loc="upper left")
ax2.set_xlabel('Month', fontsize=12)
ax2.set_ylabel('Total burned area (km^2)', fontsize=12)

#ax.hist(x, alpha=0.5, bins=100, density=True, stacked=True, label=str(cut), color=colors[i])
#ax.set_title(levels[i])
#plt.suptitle('Lithology Histogram', y=1.05, size=16)
#ax.set_xlim(50, 70); ax.set_ylim(0, 1);
plt.tight_layout();
plt.savefig(os.path.join('/qfs/projects/comexdbm/LDRD_extremes/','fire_plots','FOD&MCD64_2001_2018.png'), dpi = 300);