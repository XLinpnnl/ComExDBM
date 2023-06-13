# ComExDBM
Compound Extremes Data Benchmarking (ComExDBM)

Access by contacting xinming.lin@pnnl.gov

Project Name: Compound Extremes Data Benchmarking (ComExDBM)

Date of last update: 2023-06-12
Revision history:
2023-05-03 : Initial creation data
2023-06-12 : Cleaned up formatting
*******************************************************************************
1. Data description

Extreme weather events such as fires, heatwaves(HWs), and droughts result in significant socioeconomic
and environmental damage around the world. Mechanistic and predictive understanding of extreme weather
events are crucial for the detection, planning and response to these extremes and mitigating their 
impacts. Records of historical extreme weather events provide an important data source for understanding
present and future climate risks, but the data are sparse, unevenly distributed, and of multi-fidelity 
from multiple sources; in addition, there are many nonstandard metrics defining the levels of severity 
or impacts of extremes. In this study, we develop a benchmark data inventory of US extreme weather 
events (i.e., fires, HWs, and droughts) with daily temporal resolution and spatial resolution of 
0.5°×0.5° (~55km×55km) using data from multiple sources. 
 
The data inventory of US extreme weather events were created using heatwaves(HWs), fires, droughts, 
and meteorological variables data collected from multiple sources. The resulting datasets include 
daily temperature anomaly, heatwave labels, fire related features(i.e., fire hot spots, fire burned 
area),drought indices and co-located meteorological variables. All collected data are compiled and 
summarized to match the daily temporal resolution and 0.5°×0.5°(~55km×55km) spatial resolution. 
The output datasets are saved in netcdf files, which can be easily accessed and applied to machine 
learning (ML)-based research and would encourage ML/AI research in extreme weather and facilitate 
further work in understanding and mitigating the negative effects of these extremes.

*******************************************************************************
2. Folder structure

01_Input contains the input data information and description for all the input data from external sources 

02_Code contains all code files

03_Analysis contains the exploratory data analysis results, which include figures and docx files

04_OutputData contains seperate sub-folders for all the final output data
04_OutputData/01_HeatWave contains netcdf files for daily temperature as well as the heatwave labels 
						with spatial resolution of 0.5°×0.5°(~55km×55km)
04_OutputData/02_FireData contains netcdf files for fire related features with daily temporal resolution and  
						spatial resolution of 0.5°×0.5°(~55km×55km)
04_OutputData/03_DroughtData contains netcdf files for drought indices with daily temporal resolution and  
						spatial resolution of 0.5°×0.5°(~55km×55km)
04_OutputData/04_Meteorological_Variables contains netcdf files for meteorological variables with daily 
						temporal resolution and spatial resolution of 0.5°×0.5°(~55km×55km)
						
*******************************************************************************
3. File naming schema

File type: netcdf data files
Filename schema: ComExDBM_YYYY_DT_status.nc
Schemakey:
ComExDBM = project name
YYYY = year of data
DT = data type (i.e., HeatWave, Fires, Drought, Metvars)
status = data version 
Example filename : ComExDBM_2001_HeatWave_V01.nc

*******************************************************************************
4. File abbreviation

Project abbraviation:
ComExDBM = Compound Extremes Data Benchmarking

Filename abbraviation:
Metvars = Meteorological Variables

*******************************************************************************
5. Data sharing/access information
	a) Licenses/restrictions placed on the data:
	b) Links to the publications that cite or use the data

*******************************************************************************
6. Information about the founding source

This research has been funded by the Laboratory Directed Research and Development(LDRD) Program at 
Pacific Northwest National Laboratory (PNNL). Data processing was performed using resources available 
through Research Computing at PNNL. PNNL is operated by Battelle for the U.S. Department of Energy 
under Contract DE‐AC05‐76RL01830. 
	
