This readme file was created by by Xinming Lin. Last updated on 05/31/2023 
Access by contacting Xinming Lin at xinming.lin@pnnl.gov

Project Name: Compound Extremes Data Benchmarking (ComExDBM)
Creator : Xinming Lin; xinming.lin@pnnl.gov
Date of last update: 2023-05-31
Revision history:
2023-05-03 : Initial creation data
2023-05-31 : Cleaned up formatting

Files location : /qfs/projects/comexdbm/LDRD_extremes/ComExDBM/04_OutputData
Readme.txt location : /qfs/projects/comexdbm/LDRD_extremes/ComExDBM/04_OutputData/

*******************************************************************************
1. Folder structure, Data Type(DT) and variables description in each subfolder
	a) 01_HeatWave:  (DT - Temperature anomaly and the heatwave labels)
		This folder contains netcdf files for daily temperature as well as the heatwave labels 
		with spatial resolution of 0.5°×0.5°(~55km×55km)
		Variables included in each file:
			'date' : Date 
			'lon' : longtitude
			'lat' : latitude
			'tmax' : gridded daily maximum temperature
			'tmin' : gridded daily minimum temperature
			'tmean' : gridded daily mean temperature
			'HI02' : heat index, heatwave days are defined as a day with tmean exceeds the 95th 
					 percentile of the daily mean temperature distribution for at least 3 consecutive days
			'HI04' : heat index, heatwave days are defined as a day with tmean exceeds the 99th 
					 percentile of the daily mean temperature distribution for at least 3 consecutive days
			'HI05': heat index, heatwave days are defined as a day with tmax exceeds the 95th 
					 percentile of the daily maximum temperature distribution for at least 3 consecutive days
			'HI06': heat index, for every day in a heat wave event, the tmax must be over the 81st percentile, 
					the tmax must exceed the 97.5th percentile for at least three consecutive days within the 
					heat wave, allowing for days with tmax <97.5th percentile, the whole time period classified 
					as a heat wave must have an average tmax greater than the 97.5st percentile.
			'HI09': heat index, heatwave days are defined as a day with tmin exceeds the 95th 
					 percentile of the daily minimum temperature distribution for at least 3 consecutive days
			'HI10': heat index, for every day in a heat wave event, the tmin must be over the 81st percentile, 
					the tmin must exceed the 97.5th percentile for at least three consecutive days within the 
					heat wave, allowing for days with tmin <97.5th percentile, the whole time period classified 
					as a heat wave must have an average tmax greater than the 97.5st percentile.

	b) 02_FireData: (DT - Fire related features)
		This folder contains netcdf files for fire related features with daily temporal resolution and  
		spatial resolution of 0.5°×0.5°(~55km×55km)
		Variables included in each file:
			'date' : Date 
			'lon' : longtitude
			'lat' : latitude
			'FHS_c8c9': number of fire hotspots(i.e., number of activate fire pixels(1km×1km) within 
						each 0.5°×0.5° grid) with confidence level of 8 and 9
			'maxFRP_max_c9': daily maximum fire radiative power with the highest confidence level 
						( i.e., confidence level 9)
			'FHS_c9': number of fire hotspots(i.e., number of activate fire pixels(1km×1km) within 
					  each 0.5°×0.5° grid) the highest confidence level ( i.e., confidence level 9)
			'BA_km2': fire burned areas

	c) 03_DroughtData: (DT - Drought indices)
		This folder contains netcdf files for drought indices with daily temporal resolution and  
		spatial resolution of 0.5°×0.5°(~55km×55km)
		Variables included in each file:
			'date': Date 
			'lon': longtitude
			'lat': latitude
			'pdsi': Palmer Drought Severity Index(PDSI) 
			'pdsi_category' : the category of PDSI 
			The values of PDSI category have the following meaning:
			(
				10: 5.0 or more (extremely wet) 
				9:  4.0 to 4.99 (very wet) 
				8:  3.0 to 3.99 (moderately wet),
				7:  2.0 to 2.99 (slightly wet) 
				6:  1.0 to 1.99 (incipient wet spell) 
				5:  -0.99 to 0.99(near normal-5) 
				4:  -1.99 to -1.00 (incipient dry spell) 
				3:  -2.99 to -2.00 (mild drought) 
				2:  -3.99 to -3.00 (moderate drought) 
				1:  -4.99 to -4.00 (severe drought) 
				0:  -5.0 or less (extreme drought) 
				)
			'spi14d': the standardized precipitation index (SPI) at 14-day scale 
			'spi14d_category': category of SPI at 14-day scale 
			'spi30d' : SPI at 30-day scale 
			'spi30d_category': category of SPI at 30-day scale 
			'spi90d': SPI at 90-day scale 
			'spi90d_category': category of SPI at 90-day scale 
			'spei14d': SPEI at 14-day scale 
			'spei14d_category': category of SPEI at 14-day scale 
			'spei30d': SPEI at 30-day scale 
			'spei30d_category': category of SPEI at 30-day scale 
			'spei90d': SPEI at 90-day scale 
			'spei90d_category': category of SPEI at 90-day scale 
			The values of SPI/SPEI category have the following meaning:
			(
				10: 2.0 or more (extremely wet)
				9:  1.6 to 1.99 (very wet)
				8:  1.3 to 1.59 (moderately wet),
				7:  0.8 to 1.29 (slightly wet)
				6:  0.5 to 0.79 (incipient wet spell)
				5:  -0.49 to 0.49(near normal),
				4:  -0.79 to -0.5 (incipient dry spell)
				3:  -1.29 to -0.8 (mild drought)
				2:  -1.59 to -1.3 (moderate drought)
				1:  -1.99 to -1.6 (severe drought)
				0:  -2.0 or less (extreme drought)
				)
	d) 04_Meteorological_Variables: (DT - Meteorological variables)
		This folder contains netcdf files for meteorological variables with daily 
		temporal resolution and spatial resolution of 0.5°×0.5°(~55km×55km)
		Variables included in each file:
			'date': Date 
			'lon': longtitude
			'lat': latitude
			'air_sfc_mean', 'air_sfc_max','air_sfc_min': the daily mean/max/min of air surface temperature
			'soilm_mean', 'soilm_max','soilm_min': the daily mean/max/min of soil moisture
			'lhtfl_mean', 'lhtfl_max','lhtfl_min': the daily mean/max/min of latent heat flux
			'shtfl_mean', 'shtfl_max', 'shtfl_min': the daily mean/max/min of sensible heat flux
			'apcp_mean', 'apcp_max', 'apcp_min':the daily mean/max/min of 3-hour accumulated precipitation
			'QV_mean_250', 'QV_mean_500', 'QV_mean_850': the daily mean of specific humidity at 250/500/850 hPa
			'QV_max_250', 'QV_max_500', 'QV_max_850', : the daily max of specific humidity at 250/500/850 hPa
			'QV_min_250', 'QV_min_500', 'QV_min_850': the daily min of specific humidity at 250/500/850 hPa
			'RH_mean_250', 'RH_mean_500', 'RH_mean_850': the daily mean of relative humidity at 250/500/850 hPa
			'RH_max_250', 'RH_max_500', 'RH_max_850': the daily max of relative humidity at 250/500/850 hPa
			'RH_min_250', 'RH_min_500', 'RH_min_850' : the daily min of relative humidity at 250/500/850 hPa
			'U_mean_250', 'U_mean_250', 'U_mean_250':  the daily mean of U-wind at 250 hPa
			'U_max_500', 'U_max_500', 'U_max500':  the daily max of U-wind at 500 hPa
			'U_min_850', 'U_min_850', 'U_min_850': the daily min of U-wind at 850 hPa
			'V_mean_250', 'V_mean_250', 'V_mean_250':  the daily mean of V-wind at 250 hPa
			'V_max_500', 'V_max_500', 'V_max_500':  the daily max of V-wind at 500 hPa
			'V_min_850', 'V_min_850', 'V_min_850':  the daily min of V-wind at 850 hPa
			'BCOC_mean', 'BCOC_max', 'BCOC_min':  the daily mean/max/min of black carbon and organic carbon aerosols (BC+OC) 

*******************************************************************************
2. File naming schema

File type: netcdf data files
Filename schema: ComExDBM_YYYY_DT_status.nc
Schemakey:
ComExDBM = project name
YYYY = year of data
DT = data type (i.e., HeatWave, Fires, Drought, Metvars)
status = data version 
Example filename : ComExDBM_2001_HeatWave_V01.nc

*******************************************************************************
3. File abbreviation

Project abbraviation:
ComExDBM = Compound Extremes Data Benchmarking

Filename abbraviation:
Metvars = Meteorological Variables