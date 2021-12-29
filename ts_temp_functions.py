# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
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


#%%
#These object contains the functinos to find the nearest station. 
#This object takes as an input, a Point(lon,lat) that is the 
#reference point of interest
#The only other input, printupdate, is only used to display status updates for the script
#If false, no updates are given
#The outputs of interest are as follows:
    #self.refpoint is the reference
    #self.id_closest_st gives the ID of the nearest station
    #self.name_closest_st gives its name
    #self.miles_from_closest_st gives its distance
    #self.closest_ten gives a list of the 10 nearest stations
    
class FindStation :
    def __init__(self,point,printudpate=False):
                
        self.display=  printudpate   

        self.refpoint=point
        
        #This runs the functino and creates a variable list that is the miles distance from the ref point.
        self.closest_ten = self.nearest_station()
        
        
    #This function takes as an input, a Point(lon,lat) that is the 
    #reference point of interest
    #It then searches for the nearest point in our list of points,
    #searching distance by miles distance
    #It returns the top 10 nearest stations, ordered
    #This is returned in the foramt of a dataframe.
    #The first row of the returned dataframe is the reference point. 
    
    def nearest_station(self):
       pt=self.refpoint
       #DATAPATH is just where ghcnd-stations and ghcnd-inventory are located on your hard drive.
       stDATAPATH = "C:\\ts_big_data_files\\"
       #This initializes the actual data to search, from the list of stations
       #It only reads in the columns which are necessary to search by distance. 
       #It also searches only for the TMAX values. 
       df = pd.read_csv(stDATAPATH+'ghcnd_station_master_ts_tmax.csv')[['ID','Name','Latitude','Longitude']].drop_duplicates()
       
       #This steps strips out lat and lon values that are not nearby, reducing the number
       #of distance computations required.
       #Using this limiter reduced the time to run this script by an entire second,
       #from 1.2 seconds when searching all 40,000 rows with TMAX,
       #to 0.1 seconds, when searching within 0.25 lat /lon degrees
       bar=1
       df=df[df['Longitude']>pt.x-bar]
       df=df[df['Longitude']<pt.x+bar]
       df=df[df['Latitude']>pt.y-bar]
       df=df[df['Latitude']<pt.y+bar]
       
       #If the resulting trimmed dataframe is too small, you could introdcue errors
       #in the search. So, if the numebr is lses than a given value, 
       #then it reloads and expands the search.
       if len(df) <40:
           df = pd.read_csv(stDATAPATH+'ghcnd_station_master_ts_tmax.csv')[['ID','Name','Latitude','Longitude']].drop_duplicates()
           bar1=10
           df=df[df['Longitude']>pt.x-bar1]
           df=df[df['Longitude']<pt.x+bar1]
           df=df[df['Latitude']>pt.y-bar1]
           df=df[df['Latitude']<pt.y+bar1]
       
        #Prints th enumber of stations being searched.
       if self.display: print("Searching closest station among "+str(len(df))+" stations.")
       #This then transforms the initial data series into a GeoSeries of points   
       t1=gpd.GeoSeries(gpd.points_from_xy(df.Longitude, df.Latitude))
       #Then evaluates the distance between pt and each entry in the data series
       t1= t1.distance(pt)
       #Then sorts all distances, to find closest
       t2=t1.sort_values()
       #And selects the closest values
       #20 are selected, since in the next stage, we re-sort by actual distance (miles)
       #this is relevant because the distance, calculated above as cartesian distance using lat/lon,
       #is not accurate, and the distance in miles may very enough to change which stations are closest
       otw=df.iloc[t2.index[0:20]]
       
       #This creates a dataseries which calculates the distancef rom the ref "pt"
       #to all of the nearest 20 stations.
       d1er=[]
       for i in range(0,len(otw)):
           d1er.append(self.ts_latlon_distance([pt.y,pt.x],otw.iloc[i][['Latitude','Longitude']]))
          
       #This appends the data series to the final return dataframe
       otw['Miles_from_Ref']=d1er
       #and then sorts by miles from ref, to find the closest stations
       returner=otw.sort_values(by='Miles_from_Ref')
       
       #Then return the closest 10 stations.
       returner=returner[0:10]
       
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
       #print("Time taken with stripping extraneous rows")
       #executionTime = (time.time() - startTime)
       #print('Execution time in seconds: ' + str(executionTime))
       
       final = refdf.append(returner)
       #Prints the output, if that is selected.
       if(self.display):
           print("Here's the result.")
           print(final)
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
         
#%%
#This object contains the functions to read in key data about a station.
#Its input is a station ID. 
class StationReader:
    def __init__(self,stationID):
        self.stid = stationID
        
        #This YAML file contains a great deal of static information, 
        #such as directory information. 
        yaml_file = open("load_stats_static.yaml")
        self.yaml = yaml.load(yaml_file, Loader=yaml.FullLoader)
        
        
        self.datapath= self.yaml['DATAPATH']

        
    #This loads the station, using the station ID which is passed.
    def load_station(self,station):
        
        startTime = time.time()
        if station == None:
            station=self.stid
        #This assigns the filename where station info is located.
        filename=self.datapath+str(station)+'.dly'
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
        
        #This is to test whether the prelim, read-in step, helps or hinder the read time
        
        print("Stripping out extraneous rows before reading everything in...")
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
            
 
            
        self.station_data=inter0
        executionTime = (time.time() - startTime)
        print('Time taken to load station data using load_station: ' + str(executionTime))
        return inter0
