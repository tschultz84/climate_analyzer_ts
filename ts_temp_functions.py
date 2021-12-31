# -*- coding: utf-8 -*-

import shapely
from shapely.geometry import Point,MultiPoint,Polygon
from shapely.ops import nearest_points
import geopandas as gpd
import pandas as pd
from math import radians, cos, sin, asin, sqrt
import time
import os
import sys
import yaml
import numpy as np
from datetime import date


#%%
"""
This object has functions to identify the station nearest to a reference lat, lon.
The inputs are as follows:
*point is a list, in the form [latitude,longitude], which is the reference point
(data will be lodded for the station closest to these coordinations, unless stid is not None)
*numsts is the number of stations returned during the search, default = 10
*stid is a manually entered station. If stid is not None, then point is overriden, 
and the specific station ID is instead searched for. 
#autorun runs all the functions automatically. 
*printupdate, is only used to display status updates for the script
#If false, no updates are printed

#The outputs of interest are as follows:
    #self.refpoint is the reference
    #self.id_closest_st gives the ID of the nearest station
    #self.name_closest_st gives its name
    #self.miles_from_closest_st gives its distance
    #self.closest_stations gives a list of the  nearest stations
    self.station_data is a dataframe ahving loaded all the statin's data. '
"""    
class LoadStation :
    def __init__(self,point,numsts=10,autorun=True,printudpate=False):
                
        self.display=  printudpate   

        self.numstats = numsts

        self.refpoint=Point(point[1],point[0])
        
        #This YAML file contains a great deal of static information, 
        #such as directory information. 
        yaml_file = open("load_stats_static.yaml")
        self.yaml = yaml.load(yaml_file, Loader=yaml.FullLoader)
        
        #This automatically runs the functions.
        if autorun==True:
            self.closest_stations = self.nearest_station()
            self.station_data=self.load_station(self.closest_stations.iloc[1,0])
            self.station_data_clean=self.StationDataCleaner()
            self.station_data_clean=self.calculate_tmid(self.station_data_clean)
    """    
    #This function as a default takes the LoadStation object's
    #reference point of interest and searches for the nearest point in our list of points,
    #searching distance by miles distance
    #It returns a list of the nearest stations, ordered
    #This is returned in the foramt of a dataframe.
    #The first row of the returned dataframe is the reference point. 
    The input 'num' is the number of closest stations returned
    
    """
    def nearest_station(self):
       pt=self.refpoint
       if self.display == True : startTime = time.time()
       #DATAPATH is just where ghcnd-stations and ghcnd-inventory are located on your hard drive.
       #stDATAPATH = "C:\\ts_big_data_files\\"
       #This initializes the actual data to search, from the list of stations
       #It only reads in the columns which are necessary to search by distance. 
       #It also searches only for the TMAX values. 
       #df = pd.read_csv(stDATAPATH+'ghcnd_station_master_ts_tmax.csv')[['ID','Name','Latitude','Longitude']].drop_duplicates()
       df = pd.read_csv(self.yaml['STATIONMETA'],
                        dtype={'Firstyear': np.int64,'Lastyear': np.int64})[['ID','Name','Latitude','Longitude','Firstyear','Lastyear']].drop_duplicates()
       
       #These lines strip out stations where there is no recent data (from within the last year), 
       recentyear=date.today().year-1
       df = df[df['Lastyear']>=recentyear] 
       #and then strips out stations for which data is only very very recent. 
       baseyear=self.yaml['BASEYEAR']
       df = df[df['Firstyear']<=baseyear] 
       
       #Now drops the year info, since we don't need it.
       df=df[['ID','Name','Latitude','Longitude']]
           
       
       #This steps strips out lat and lon values that are not nearby, reducing the number
       #of distance computations required.
       #Using this limiter reduced the time to run this script by an entire second,
       #from 1.2 seconds when searching all 40,000 rows with TMAX,
       #to 0.1 seconds, when searching within 0.25 lat /lon degrees
       bar=self.yaml['SEARCH_RADIUS']
       df=df[df['Longitude']>pt.x-bar]
       df=df[df['Longitude']<pt.x+bar]
       df=df[df['Latitude']>pt.y-bar]
       df=df[df['Latitude']<pt.y+bar]
       
       #If the resulting trimmed dataframe is too small, you could introdcue errors
       #in the search. So, if the numebr is lses than a given value, 
       #then it reloads and expands the search.
       if len(df) < self.numstats*4:
           df = pd.read_csv(self.yaml['STATIONMETA'])[['ID','Name','Latitude','Longitude']].drop_duplicates()
           bar1=10
           df=df[df['Longitude']>pt.x-bar1]
           df=df[df['Longitude']<pt.x+bar1]
           df=df[df['Latitude']>pt.y-bar1]
           df=df[df['Latitude']<pt.y+bar1]
       
        #Prints th enumber of stations being searched.
       if self.display: 
           print("Searching closest station among "
                              +str(len(df))+" stations within "+str(self.yaml['SEARCH_RADIUS'])+" degrees of the reference.")
           print("This only includes stations with data available more recently than "+str(recentyear)+" and before "+str(self.yaml['BASEYEAR']))
       #This then transforms the initial data series into a GeoSeries of points   
       t1=gpd.GeoSeries(gpd.points_from_xy(df.Longitude, df.Latitude))
       #Then evaluates the distance between pt and each entry in the data series
       t1= t1.distance(pt)
       #Then sorts all distances, to find closest
       t2=t1.sort_values()
       #And selects the closest values
       #4*num are selected, since in the next stage, we re-sort by actual distance (miles)
       #this is relevant because the distance, calculated above as cartesian distance using lat/lon,
       #is not accurate, and the distance in miles may very enough to change which stations are closest
       otw=df.iloc[t2.index[0:4*self.numstats]]
       
       #This creates a dataseries which calculates the distancef rom the ref "pt"
       #to all of the nearest 20 stations.
       d1er=[]
       for i in range(0,len(otw)):
           d1er.append(self.ts_latlon_distance([pt.y,pt.x],otw.iloc[i][['Latitude','Longitude']]))
          
       #This appends the data series to the final return dataframe
       otw['Miles_from_Ref']=d1er
       #and then sorts by miles from ref, to find the closest stations
       returner=otw.sort_values(by='Miles_from_Ref')
       
       #Then return the closest numstats stations.
       returner=returner[0:self.numstats]
       
       #This sets a variable which is the ID Of the nearest station
       self.id_closest_st = returner.iloc[0,0]
       #Then a separate variable that is the name. 
       self.name_closest_st = returner.iloc[0,1]
       #Then a separate variable that is the miles distance from the ref point. 
       self.miles_from_closest_st = returner.iloc[0,4]
       
       #This creates the frist row of the return dataframe.
       refdf= pd.DataFrame(
           {'ID':['REFPT'],
            'Name':["Reference Location"],
            'Latitude':[pt.y],
            'Longitude':[pt.x],
            'Miles_from_Ref':[0]}
           )
       
       

       final = refdf.append(returner)
       #Prints the output, if that is selected.
       if(self.display):
           executionTime = (time.time() - startTime)
           print("Time taken to locate nearest "+str(self.numstats)+" weather stations." + str(executionTime))
    
       return final
    

    #This function returns the miles between two poitns
    #def ts_latlon_distance(lat1, lat2, lon1, lon2):
    def ts_latlon_distance(self,latlon1,latlon2):
        # radians which converts from degrees to radians.
        lon1 = radians(latlon1[1])
        lon2 = radians(latlon2[1])
        lat1 = radians(latlon1[0])
        lat2 = radians(latlon2[0])
        
        # Haversine formula
        dlon = lon2 - lon1
        dlat = lat2 - lat1
        a = sin(dlat / 2)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2)**2
     
        c = 2 * asin(sqrt(a))
        
        # Radius of earth in kilometers. Use 3956 for miles
        r = 3956
          
        # calculate the result, returned in miles
        return(c * r)
         
    #This loads the station.
    #As default, station = None, and the id_closest_st for thsi object is run.
    #The station value must be a string. 
    def load_station(self,station=None):
        #This timer is used to check how long it takes to run the station read in function. 
        if self.display == True :
           startTime = time.time()
        #If statino is none, then defaults to the station closest to the ref point. 
        if station == None:
            station=self.id_closest_st
        
        #This assigns the filename where station info is located.
        filename=self.yaml['DATAPATH']+str(station)+'.dly'
        #This checks the file exists and breask if it does not.
        #(If it doesn't exist, the submitted station ID was in error.)
        if os.path.exists(filename) == False:
            print("Bad Station ID. The file called "+filename+" does not exist.")
            print("Please check your station ID "+str(station)+" and re-submit.")
            sys.exit("Break Error in load_station of StationReader: Bad Station ID.")
        #You have designed this script so it only reads in the lines where data is useful
        #and excludes other elements. 
        #When testing reading in the very largest data file, you would typically save
        #0.02 to 0.06 seconds during loading. Worth every bit...
        
        #This defines the pre-loading data, i.e. which columns are used to skip rows.
        shortcols = [ (17, 21)]
        shortnames =  np.array([ 'Element'])
        
        #This then reads the file, but only the columns which generate the rows to skip. 
        prelim = pd.DataFrame()    
        prelim=pd.read_fwf(filename,colspecs=shortcols)
        prelim.columns=shortnames
        #This generates the rows to skip.   
        #Currently you are keeping only TMAX and TMIN
        tokeep=['TMAX','TMIN']
        skiprows1 = prelim[(~prelim['Element'].isin(tokeep))].index
        
        
        #Then, the data is read in.
        inter0 = pd.DataFrame()      
        #This then Reads out the next file.  
        inter0 = pd.read_fwf(filename,colspecs=self.yaml['datacolnums'],skiprows=skiprows1+1)
        inter0.columns=self.yaml['datacolnames']
            
    
        #This is a dataframe of all the data.
        self.station_data=inter0
        
        if self.display == True :
            executionTime = (time.time() - startTime)
            print('Time taken to load data for +'+str(station)+' using load_station: ' + str(executionTime))
        return inter0
 
    
    """
      This takes in data for single weather station - in the format of output from
      load_station
      and cleans it up. 
      The only input is stationd, which, if =None, just defaults to using data frmo
      self.station_data
      Currently these cleaning operations are completed:
          *Change -9999 values to np.nan
    """
    def StationDataCleaner(self,stationd=None):
        #If station=d is none, just default to using self.station_data
        if stationd==None: stationd = self.station_data
        #Does the swap of -9999 for np.nan
        returner = stationd.replace(-9999,np.nan) 
        
        #Converting frmo tenths of degrees C to C.
        returner.iloc[:,-31:] = returner.iloc[:,-31:]/10
        
        #Convert from Celsius to Farnehit
        returner.iloc[:,-31:] = (returner.iloc[:,-31:]*9/5)+32
        
        return returner
    
    def calculate_tmid(self,stationd):
        if self.display==True: startTime = time.time()
        
        #Finds the earliest and latest years.
        yrmin = stationd['Year'].min()
        yrmax = stationd['Year'].max()


        stationd=stationd.sort_values(by=['Year','Month'])
        statID=stationd.iloc[0,0]
        tmiddata=pd.DataFrame()

        #The ' k' loop is over all elements.
        for k in np.arange(0,len(stationd)-1):
       # for k in np.arange(0,10):
            #This initializes the temporary dataframes. 
            calc=pd.DataFrame()
            calcint=pd.DataFrame()
            #This extracts the data for the right year and months.
            calcmax = stationd.iloc[k]
            calcmin = stationd.iloc[k+1]
            
            if calcmax['Element'] == "TMAX" and calcmin['Element'] == "TMIN":
                if calcmax['Year'] == calcmin['Year'] and calcmax['Month'] == calcmin['Month']:
                      
                    #This combines TMAX and TMIN values, such that the average can be taken. 
                    calcint1=pd.concat([calcmax,calcmin],axis=1)
                   # returner=calcint1
                    #This calcualtes the TMID values, which are the average of 
                    #TMAX and TMIN.
                    #It is importnat that skipna=FAlse, which ensures that
                    #if either tMAX or TMIN is NaN, then so is the result for TMID.
                    calcint = calcint1.iloc[4:].sum(axis=1,skipna=False)/2
                    #print(calcint)
                    #This inserts back in the missing collumns. 
                    calcint['Year']=calcmax['Year']
                    calcint['Month']=calcmax['Month']
                    calcint['Element']='TMID'
                    calcint['Station_ID']=statID
                    #And this appends the new row to the returned dataframe. 
                    tmiddata=tmiddata.append(calcint,ignore_index=True)
        returner=stationd.append(tmiddata)

        if self.display==True: 
            executionTime = (time.time() - startTime)
            print("Time taken to create TMID by going through rows " + str(executionTime))
        return returner
        
#All analysis fucntions for a station are done in this object.
#The input is stationdata, which shoudl be self.station_data_claned from the StationLoad object
#the refperiod is a 2 element list of start and end, reference period of analysis
#autorun simply immediately runs the KPI key_metrics function if True 

class StationAnalyzer :
    def __init__(self,stationdata,refperiod=[np.datetime64('2020-01-31'),np.datetime64('2020-12-31')],autorun=True):
        
        #separates out the different variables ofinterest.
        self.station_data_tmid=stationdata[stationdata['Element']=='TMID']
        self.station_data_tmax=stationdata[stationdata['Element']=='TMAX']
        self.station_data_tmin=stationdata[stationdata['Element']=='TMIN']
        
        self.refstart=refperiod[0]
        self.refend=refperiod[1]
        
        #Runs the key_metrics function if desired
        if autorun==True:
            self.kpi=self.key_metrics()
            #print(self.kpi)
        
    #This creates a table (pd dataframe) of the key metrics of interest
    def key_metrics(self):
        
        whole_tmid_mean=np.nanmean(self.station_data_tmid.iloc[:,-31:])
        whole_tmax=np.nanmax(self.station_data_tmax.iloc[:,-31:])
        whole_tmin=np.nanmin(self.station_data_tmin.iloc[:,-31:])
        
        #This finds the index location of the max temperature.
        maxloc=np.where(self.station_data_tmax.iloc[:,-31:]==whole_tmax)
        maxyear=  int( self.station_data_tmax.iloc[int(maxloc[0]),1]     )
        maxmonth  =  int( self.station_data_tmax.iloc[int(maxloc[0]),2]  )
        maxday  =   self.station_data_tmax.columns[int(maxloc[1])]  
        maxdate=str(maxyear)+'-'+str(maxmonth)+"-"+str(maxday[-2:])

        returner = pd.DataFrame(
            {"TMID_av_F_entire_record":[whole_tmid_mean],
             "TMAX_F_entire_record":[whole_tmax],
             "TMAX_date":[maxdate],
             "TMIN_F_entire_record":[whole_tmin]
             
             }
            )
        return returner           