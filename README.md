# Package Author and Use
This package was originally written by Tobias Schultz in 2022. First relase on May 10, 2022.

Feel free to use as you see fit, but do attribute it to me if you ever publish or write anything. Thanks!

# README
This file describes key environment variables, and the needed program
inputs and outputs.

# OVERVIEW
This Python module allows you to search for the land-based weather station nearest to you. You provide the programs with a latitude and longitude,
and a time period of interest; it automatically finds the closest weather station. In addition, it filters the weather stations to find ones meeting
certain quality criteria, which you specify as inputs to the class. The output is a numpy array containing daily temperature data, meant to be analyzed by other packages. 

# About GHCND
This data is all sourced from the NOAA Dataset called Global Historical Climatology Network daily (GHCNd). To quote NOAA:
"The Global Historical Climatology Network daily (GHCNd) is an integrated database of daily climate summaries from land surface stations across the globe. 
GHCNd is made up of daily climate records from numerous sources that have been integrated and subjected to a common suite of quality assurance reviews.
GHCNd contains records from more than 100,000 stations in 180 countries and territories."

NOAA provides GHCNd data through the National Centers for Environmental Information (NCEI). See 
https://www.ncei.noaa.gov/products/land-based-station/global-historical-climatology-network-daily

# ENVIRONMENT
It should work with Python 3.9.7 and above. 

Make sure the following files are all in your root directory:
* **Station_analyzer_ts.py** contains the analysis module. 
* **Station_Loader_ts.py** contains the module to load in station data. *This is a vital module and I encourage you to read its accompanying README file - README_StationLoader.md -- to learn more about it. 
* Files in **/data** directory are vital to the correct operation of programs, especially **ghcnd_station_mast_ts_tmax.csv**. StationLoader will not work if
**data/ghcnd_station_mast_ts_tmax.csv** does not exist.
* **examples.ipynb** shows you how to run the Python modules. You don't technically need it, but I would recommend taking a look and then making a copy of this
to run other weather station analyses.
* **create_station_master_1.pynb** was used only to generate ghcnd_station_mast_ts_tmax.csv from underlying metadata in the GHCND database. You can ignore it. 

You will also need to ensure your Python programs know where to look, either by setting the current working directory using os.cwd() to the root directory or adding a new .pth file in your site-packages
directory that points there. 

You also need to ensure these packages are available in your environment: *pandas, math, time, sys, numpy, datetime, requests, warnings*

# StationAnalyzer Class
This is the vital class that is used in all analysis.

StationAnalyzer takes in this input fields:
1. **point** in the form [float,float] which is the latitude and longitude of interest to you. 
2. **printupdate**, which is TRUE or FALSE, and just indicates whether detailed processing updates are returned to your Python console. 
3. **Min_days_per_mo** is the minimum number of days of data which must be present for every month of a year for a dataset for the year to be included.  If not, then the entire year is excluded from dataset.
4. **Search_radius** is the lat/lon distance from **point** which is searched for weather stations. This is used to streamline how fast the "closest station" function works. Stations are only reviewed for quality if they are within +/- **search_radius** of **point**. 
5. **Firstyear** is the earliest year which must be present in the dataset record for the weather station
        to be loaded. If the weather station does not have some data present going back to this year, then it is rejected and another weather station will be reviewed instead.
6. **Lastbaseyear** is the last year in which the baseline will be calculated. The baseline period is defined as the period between **firstyear** and **lastbaseyear**. This variable is important for now because there must be **basenoyears** of data available in the baseline period. 
7. **Basenoyears** is the number of years required in the baseline period. For example, if **basenoyears == 30**, then there must be 30 complete years of data between **firstyear** and **lastbaseyear**.
8. **Min_recent_years** is the minimum number of years before the present which must be present for the statin to be loaded.
9. **Required_trend_years** is the minimum number of years in the last 30 required to calculate a trend.   Weather stations without this many years present in the last 30 years of data will be discarded.

The function called **change_ref_dates** is also critical, as you use it to change your reference period information. This function takes in as its first argument *year_series*, a list of [year_start,year_end] which is the beginning and end of the period. The 2nd argument is *date_series*, list of [[month, day],[month, day]] which is the beginning and end of a subannual reference. If *dtype(date_series)* != list, then it assumes the beginning and end are the whole year. 

Two other key functions are **key_metrics**, which generatesthe tables called **key_data** and **key_stats**, and **key_charts**, which outputs charts showing histograms of the daily temperatures, trend in temperature over the entire historical period (highlighting the reference, baseline, and trends) and monthyl average temperatures over the entire record. 

# Attributes of StationAnalyzer
StationAnalyzer is a class. After you run it, the following object attributes are available. 
1) **StationAnalyzer.stationobj** is an object which contains all of the attributes associated with the StationLoader class, which includes:
   
    a) **StationAnalyzer.stationobj.station_data** is a Nx7 numpy array that is the actual temperature data for the selected weather station. Each row is a day; and it has these 7 columns: *[Year, Month, Day, Day of Year, TMAX (max. temperature), TMIN (min. temperature), TMID (midpoint temperature between TMAX and TMIN)]*.

    b) **StationAnalyzer.stationobj.station_information** is a pandas DataFrame containing meta information for the selected weather station, like the station name, ID, lat/lon, etc.  '
    c) **StationAnalyzer.stationobj.closest_stations** is a pandas DataFrame which lists the weather stations closest to **point** (closest based on number of miles between) values, along with their index numbers, names, locations, and distances from point.

    d) **StationAnalyzer.stationobj.station_filters** is a pandas DataFrame which summarizes all the filters you entered as inputs (e.g., firstyear, search_radius, etc).
2) **StationAnalyzer.all_data_long** and **StationAnalyzer.tmid_array**, which are respectively [Nx7] and [Nx5] numpy arrays. They contain these columns: [Year,Month,Day,Day of Year,TMAX,TMIN,TMID,and Julian Date], although **StationAnalyzer.tmid_array** excludes TMID. Every row is a day.
3) **StationAnalyzer.period_information** is a pandas Data series including information on the period analyzed (baseline period, reference period, etc.)
4) **StationAnalyzer.all_years_mean** and **StationAnalyzer.all_months_mean** are respectively averages across the selected subset of time.
5) **StationAnalyzer.key_data** and **StationAnalyzer.key_stats** respectively tables containing a summary of the climatological information, and temperature changes between the reference and baseline periods and in the most recent 30 years. 
