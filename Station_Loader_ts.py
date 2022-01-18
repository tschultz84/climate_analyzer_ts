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
import requests


#%%
"""
This object has functions to identify the station nearest to a reference lat, lon.
The inputs are as follows:
*point is a list, in the form [latitude,longitude], which is the reference point
(data will be lodded for the station closest to these coordinations, unless stid is not None)

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
    self.station_data is a dataframe ahving loaded all the statin's data, a numpy array 
     with these columns: [Year, Month, Day, Day of Year, TMAX, TMIN, TMID]
"""    
class LoadStation :
    def __init__(self,point,printudpate=False):
                
        self.display=  printudpate   
        #If a point (lat, lon) is passed, then create a point. 
        if (type(point) != str) and (len(point)==2):
            self.refpoint=point
        
        #This YAML file contains a great deal of static information, 
        #such as directory information. 
        yaml_dir="C:\\Users\\14154\\OneDrive\\Python\\climate_mapper\\python\\climate_analyzer_ts\\"
        yaml_file = open(yaml_dir+"load_stats_static.yaml")
        self.yaml = yaml.load(yaml_file, Loader=yaml.FullLoader)
        """
        The following few lines automatically run all the functions at initiation.
        """
        if self.display: startTime = time.time() 
        #If point is a string, you skip finding any station, and go righ o the loader. 
        if type(point) == str:
            self.station_data=self.load_station(point)
        
        #If point is not a string, but a point, run the function to find
        #the nearest station.     
        if type(point) != str:    
            self.closest_stations = self.nearest_station(point)
            if(self.display): 
                print("Time to locate nearest weather stations:" + str(time.time() - startTime))
                startTime = time.time() 
            self.station_data=self.load_station(self.id_closest_st)
        if(self.display): 
            print("Time to load this station data:" + str(time.time() - startTime))
            startTime = time.time() 
        
        self.station_data=self.StationDataCleaner(self.station_data)
        
        if(self.display): 
            print("Time to clean up this data:" + str(time.time() - startTime))
    """    
    #This function as a default takes the LoadStation object's
    #reference point of interest and searches for the nearest point in our list of points,
    #searching distance by miles distance
    #It returns a list of the nearest stations, ordered
    #This is returned in the foramt of a dataframe.
    #The first row of the returned dataframe is the reference point. 
    The input 'num' is the number of closest stations returned
    
    """
    def nearest_station(self,point):
       pt=point
       
       #DATAPATH is just where ghcnd-stations and ghcnd-inventory are located on your hard drive.
       #stDATAPATH = "C:\\ts_big_data_files\\"
       #This initializes the actual data to search, from the list of stations
       #It only reads in the columns which are necessary to search by distance. 
       #It also searches only for the TMAX values. 
       df1 = pd.read_csv(self.yaml['STATIONMETA'],
                        dtype={'Firstyear': np.int64,'Lastyear': np.int64})[['ID','Name','Latitude','Longitude','Firstyear','Lastyear']]
       
       #These lines strip out stations where there is no recent data (from within the last year), 
       recentyear=date.today().year-1
       df = df1[df1['Lastyear']>=recentyear] 
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
       df=df[df['Longitude']>point[1]-bar]
       df=df[df['Longitude']<point[1]+bar]
       df=df[df['Latitude']>point[0]-bar]
       df=df[df['Latitude']<point[0]+bar]

        #Prints th enumber of stations being searched.
       if self.display: 
           print("Searching closest station among "
                              +str(len(df))+" stations within "+str(self.yaml['SEARCH_RADIUS'])+" degrees of the reference.")
           print("This only includes stations with data available more recently than "+str(recentyear)+" and before "+str(self.yaml['BASEYEAR']))
       """
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
       otw=df.iloc[t2.index[0:40]]
       """
       #This creates a dataseries which calculates the distancef rom the ref "pt"
       #to all of the nearest 20 stations.
       d1er=[]
       for i in range(0,len(df)):
           d1er.append(self.ts_latlon_distance([point[0],point[1]],df.iloc[i][['Latitude','Longitude']]))
       
       #This appends the data series to the final return dataframe
       df['Miles_from_Ref']=d1er
       
       #and then sorts by miles from ref, to find the closest stations       
       returner=df.sort_values(by='Miles_from_Ref')
       
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
            'Latitude':[point[0]],
            'Longitude':[point[1]],
            'Miles_from_Ref':[0]}
           )
       
       

       final = refdf.append(returner)
       #Prints the output, if that is selected.
       
    
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

        #If statino is none, then defaults to the station closest to the ref point. 
        if station == None:
            station=self.id_closest_st
        """
        #This assigns the filename where station info is located.
        filename=self.yaml['DATAPATH']+str(station)+'.dly'
        
        #This checks the file exists and breask if it does not.
        #(If it doesn't exist, the submitted station ID was in error.)
        if os.path.exists(filename) == False:
            print("Bad Station ID. The file called "+filename+" does not exist.")
            print("Please check your station ID "+str(station)+" and re-submit.")
            sys.exit("Break Error in load_station of StationReader: Bad Station ID.")
        """
        #Downloading the CSV file straight rom NOAA.
        #First, assigning the URL Name.
        csv_url = self.yaml['NOAA_URL']+str(station)+'.csv'
        
        if self.display == True: 
            #This gets the filesize.
            info=requests.head(csv_url)
            filesize = int(info.headers['Content-Length'])/1000000
            print("Beginning to download "+csv_url+' ('+str(filesize)+"MB). Please wait while the file is transferred.")
        #Then, checking i it exists. the script breaks if the file does not exist.
        response = requests.get(csv_url)
        if response.status_code!= 200:
            
            print("Bad Station ID. The file called "+csv_url+" does not exist.")
            print("Please check your station ID "+str(station)+" and re-submit.")
            sys.exit("Break Error in load_station of StationReader: Bad Station ID.")
        #Then, downloads the file. 
        inter0 = pd.read_csv(csv_url)
        """
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
        tokeep=['TMAX','TMIN','TAVG']
        skiprows1 = prelim[(~prelim['Element'].isin(tokeep))].index
        
        
        #Then, the data is read in.
        inter0 = pd.DataFrame()      
        #This then Reads out the next file.  
        inter0 = pd.read_fwf(filename,colspecs=self.yaml['datacolnums'],skiprows=skiprows1+1)
        """

        #This is a dataframe of all the data.
        #self.station_data=inter0
        
        return inter0
 
    #This function finds the Day of year, assuming every month has 31 days
    #taking in the day of month (dom) and month.
    def find_doy(self,dom,month):
        doy = dom+(month-1)*31
        return doy
    """
      This takes in data for single weather station - in the format of output from
      load_station
      and cleans it up. 
      The only input is stationd, which must be in the format of self.station_data
     It returns a numpy array version, where each row contains:
         Year, Month, Day, Day of Year, TMAX, TMIN, TMID
    """
    def StationDataCleaner(self,stationd):
       
        #Does the swap of -9999 for np.nan (I am not sure this is relevant for CSV)
        returner = stationd.replace(-9999,np.nan) 

        
        #CHeck if TAVG exists. if not, calculate it.
        if ("TAVG" in returner.columns)!=True: 
           if self.display : print("TAVG not present, I will create it.")
           returner['TMID'] = (returner['TMAX']+returner['TMIN'])/2
        #If it does exist, rename it to TMID, beacause TAVG is a milseading name.
        if ("TAVG" in returner.columns)==True: 
           if self.display : print("TAVG present, I will rename it to TMID.")
           returner = returner.rename(columns ={'TAVG':"TMID"})

        #Now, drop extraneous columns.
        returner = returner[self.yaml['KEEP_COLS']]

                    
        #Converting frmo tenths of degrees C to C.
        threetemps=["TMAX",'TMIN',"TMID"]
        returner[threetemps] = returner[threetemps]/10
        
        #Convert from Celsius to Farnehit
        returner[threetemps] = (returner[threetemps]*9/5)+32
        
        #Create dates columns appropriately.
        returner['Year']=returner['DATE'].str[:4]
        returner['Month']=returner['DATE'].str[5:7]
        returner['Day']=returner['DATE'].str[8:10]
        
        returner['Year']=pd.to_numeric(returner.Year)
        returner['Month']=pd.to_numeric(returner.Month)
        returner['Day']=pd.to_numeric(returner.Day)
        
        returner['DOY'] = self.find_doy(returner['Day'],returner['Month'])
        
        returner = returner.drop(['DATE'],axis=1)
        
        #Re-order them.
        returner1=returner.reindex(columns=['Year','Month','Day','DOY','TMAX','TMIN',"TMID"])
        #Change to numpy.
        returner1 = np.asarray(returner1)
        return returner1
        """
        #Sort values. 
        returner=returner.sort_values(by=['Year','Month'])
        
        #TURN 31 COLUMN ENTRIES, ONE FOR EVERY DAY OF YEAR, TO ONE LONG MATRIX
        #You may consider moving this into the station_data object. 
        #This allows you to number each day of the year.
        #Making it much easier to serach by days. 
        all_data=np.asarray(returner.drop(['Station_ID'],axis=1))
        self.all_data_1=all_data
        #return returner
        #You need to drop all years without 12 months of data.
        #That screws up your calculations
        #and will also bias results. 
        
        inter1=np.empty((len(all_data)*31,6))
        #inter1=np.asarray(inter1,dtype=dtyper)
        tot=0
    
        #Loop over entire length of the dataset. 
        for k in np.arange(0,len(all_data)):
            #Pull out the year, month, and then monthdata, for every row.
            year=all_data[k,0]
            month=all_data[k,1]
            monthdata=np.transpose(all_data[k,3:])
            #Turning the elemnt value into a number, which is easier to work with
            #in numpy (and faster)
            #The integer values are arbitrary. 
            ele = all_data[k,2]
            eler=-1
            if(ele == "TMIN"):
                eler=self.yaml['TMIN_INDEX']
            if(ele == "TMAX"):
                eler=self.yaml['TMAX_INDEX']
           # if(ele == "TMID"):
           #     eler=2

            #Then loop over each month, creating a enw row for every single element. 
            for i in np.arange(1,32):
                #Assign the Day of Year, which goes to more than 365 (31 * 12)
                doy = i + (month-1)*31
                
                #Then assigning each column value. 
                inter1[tot,0]=year
                inter1[tot,1]=month
                inter1[tot,2]=i
                inter1[tot,3]=doy
                inter1[tot,5]=eler
                inter1[tot,4]=monthdata[i-1]  
                tot=tot+1
        
        #This then calculates TMEAN values. 
        returner = self.calculate_tmid(inter1)
        self.all_data_np=returner
        return returner
    
    #Calculates TMID. 
    #Station d must be in the format of self.all_data_np. 
    def calculate_tmid(self,stationd):
        if self.display==True: startTime = time.time()
        
        #Breaks apart the TMAX and TMIN values. 
        tmax=stationd[np.where(stationd[:,5]==self.yaml["TMAX_INDEX"])]
        tmin=stationd[np.where(stationd[:,5]==self.yaml["TMIN_INDEX"])]

        #THen loops over every entry, finds the mtaching dates int he other array,
        #then, if they are teh same, adds it to the returned array.
        
        outer= []
        for i in np.arange(0,len(tmax)):
            year=tmax[i,0]
            doy = tmax[i,3]
            index = np.where((tmin[:,0]==year)&(tmin[:,3]==doy))
            if np.shape(index)[1]>0:
                tminval = tmin[index,4][0][0]
                outer.append([tmax[i,0],tmax[i,1],tmax[i,2],tmax[i,3],tmax[i,4],tminval])
        
        outer=np.asarray(outer)
        #Evaluates the means. 
        means = np.nanmean(np.array([outer[:,4],outer[:,5]]),axis=0)
        #Grabs the dates. 
        dates=outer[:,0:4]
        #Then combines verything into a single matrix to be returned. 
        returner1 = np.insert(dates,4,means,axis=1)
        returner1=np.insert(returner1,5,self.yaml['TMID_INDEX'],axis=1) 
        
        #These steps then scrub out any years in TMID day which
        #do not include 12 months of data
        #This is improtant, since if there are not all moths included
        #It would bias results like yearly averages.
        #First, assign a new variable, tmidex, to be ecleaned up.
        tmidex=returner1
        
        wheres = []
        #This loops over all of the unique years in the dataset.
        uniqueyears=np.unique(tmidex[:,0])
        for i in uniqueyears:
            #First, finds the index values where the year is equal to year i.
            index=np.where(tmidex[:,0]==i)
            #Then, finds the unique months for this year, which shoudl be1 through 12.
            months=np.unique(tmidex[index][:,1])
            if len(months) ==12:
                    #Only if the year has 12 months, is the year included in the final dataset. 
                 wheres.append(index)
        wheres = np.asarray(wheres).flatten()
        tmidex1 = tmidex[wheres]
        
        #This ten combines the TMID, TMAX, and TMIN, into one big matrix. 
        returner1=np.append(tmidex,tmax,axis=0)
        #returner1=np.append(tmidex1,tmax,axis=0)
        returner1=np.append(returner1,tmin,axis=0)
        if self.display==True: 
            executionTime = (time.time() - startTime)
            print("Time taken to create TMID using NUMPY " + str(executionTime))
        
        return returner1
    """
              