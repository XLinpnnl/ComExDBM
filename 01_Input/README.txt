This readme file was created by by Xinming Lin. Last updated on 05/31/2023 
Access by contacting Xinming Lin at xinming.lin@pnnl.gov

Project Name: Compound Extremes Data Benchmarking (ComExDBM)
Creator : Xinming Lin; xinming.lin@pnnl.gov
Date of last update: 2023-05-31
Revision history:
2023-05-03 : Initial creation data
2023-05-31 : Cleaned up formatting

Input files location : /qfs/projects/comexdbm/LDRD_extremes/ComExDBM/01_Input
Readme.txt location : /qfs/projects/comexdbm/LDRD_extremes/ComExDBM/01_Input/

*******************************************************************************
Input Data Introduction:

The ComExDBM data product is produced using fires, heatwaves(HWs), droughts and meteorological 
data collected through different sources. Here is the detailed description for each data source:

1)  Fire data from the Terra Moderate Resolution Imaging Spectroradiometer (MODIS) 
   
    a) MODIS Thermal Anomalies and Fire Daily (MOD14A1) Version 6 
	The datasets are generated at ~1 kilometer (km) spatial resolution and daily temporal resolution. 
	The variables include the fire mask, pixel quality indicators, maximum fire radiative power (MaxFRP), 
	and the position of the fire pixel within the scan. Individual 1-km pixels (grids) are assigned to 
	one of nine fire mask pixel classes, which indicate the different confidence levels of fire occurrence. 
	In this study, we only use the fire pixels (grids) with the highest confidence level of fire occurrence 
	to summarize the daily fire features which include mean and maximum of MaxFRP, total number of fire hotspots (FHS). 
   
    MOD14A1 can be downloaded at https://lpdaac.usgs.gov/products/mod14a1v006/
   
    b) MODIS Burned Area Product (MCD64A1) Version 6
	MCD64A1 is a burned area product generated at 500m spatial resolution and monthly temporal resolution. 
	The Fire Events Delineation (FIRED), an event-delineation algorithm, has been used to derive the daily fire 
	burned area from the MODIS MCD64 burned area product for the coterminous US (CONUS) from 2001 to 2020.
   
	MCD64A1 can be downloaded at https://lpdaac.usgs.gov/products/mcd64a1v006/
	FIREDpy: Python Command Line Interface for classifying fire events from MODIS MCD64A1 Version 6
			https://github.com/earthlab/firedpy
   
2)  CPC Global Unified Temperature
	The gauge-based Climate Prediction Center (CPC) Global Unified daily gridded temperature data are 
	provided by NOAA Physical Sciences Laboratory.The data is global telecommunications system (GTS) data 
	and is gridded using the Shepard Algorithm. The spatial resolution of this data is 0.5°×0.5°(~55km×55km)
	All HW events in this study are summarized using the CPC Global Unified daily temperature anomaly data.
	
	CPC Global Unified Temperature can be downloaded at https://psl.noaa.gov/data/gridded/data.cpc.globaltemp.html
	
3)  Drought data from Gridded Surface Meteorological (gridMET)
    The drought indices are obtained from Gridded Surface Meteorological (gridMET) dataset with 4-km spatial resolution
	and 5-day temporal resolution. The drought indices provided include:
		standardized precipitation index (SPI), 
		standardized precipitation evapotranspiration index (SPEI), 
		Palmer Drought Severity Index (PDSI) etc. 
	The SPI and SPEI indices are supplied on different time scales corresponding to the time aggregation of precipitation, 
	reference evapotranspiration, and precipitation minus reference evapotranspiration, respectively. 
	
	The drought indices gridMET can be downloaded at https://www.climatologylab.org/gridmet.html
	
4)  Meteorological variables
	a) Meteorological variables from the North American Regional Reanalysis (NARR)
	Meteorological variables from NARR is available every 3-hour at 32-km horizontal grid spacing and 45 vertical layers
	Variables from NARR are showned as follows:
		precipitation (P3h) , soil moisture (SM), latent heat flux (LHF), and sensible heat flux (SH) at 250 and 850 hPa
	
	The NARR data can be downloaded at https://psl.noaa.gov/data/gridded/data.narr.html
	
	b) Meteorological variables from the Modern-Era Retrospective Analysis for Research and Applications Version 2 (MERRA-2)
	Meteorological variables from MERRA-2 is available every 3-hour at an approximate spatial resolution of 0.5° × 0.625° 
	and 72 hybrid-eta levels. 
	Variables from MERRA-2 are showned as follows:
		specific humidity (Qv) at 250 and 850 hPa, 
		relative humidity (RH), U-wind and V-wind at 250, 500, and 850 hPa, 
		carbonaceous aerosols, i.e., black carbon and organic carbon aerosols (BC+OC)
	
	The MERRA-2 data can be downloaded at https://disc.gsfc.nasa.gov/datasets/M2I3NPASM_5.12.4/summary
 
*******************************************************************************

Reference:

Giglio, L., C. Justice. MOD14A1 MODIS/Terra Thermal Anomalies/Fire Daily L3 Global 
1km SIN Grid V006. 2015, distributed by NASA EOSDIS Land Processes DAAC, 
https://doi.org/10.5067/MODIS/MOD14A1.006.

Giglio, L., C. Justice, L. Boschetti, D. Roy. MCD64A1 MODIS/Terra+Aqua Burned Area 
Monthly L3 Global 500m SIN Grid V006. 2015, distributed by NASA EOSDIS Land Processes 
DAAC, https://doi.org/10.5067/MODIS/MCD64A1.006

Balch, Jennifer K., Lise A. St. Denis, Adam L. Mahood, Nathan P. Mietkiewicz, 
Travis M. Williams, Joe McGlinchy, and Maxwell C. Cook. "FIRED (Fire Events Delineation):
an open, flexible algorithm and database of US fire events derived from the MODIS burned 
area product (2001–2019)." Remote Sensing 12, no. 21 (2020): 3498. 
https://doi.org/10.3390/rs12213498