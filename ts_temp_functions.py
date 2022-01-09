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
import matplotlib.pyplot as plt


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
      #      self.station_data_clean=self.calculate_tmid(self.station_data_clean)
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
       df=df[df['Longitude']>pt.x-bar]
       df=df[df['Longitude']<pt.x+bar]
       df=df[df['Latitude']>pt.y-bar]
       df=df[df['Latitude']<pt.y+bar]

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
        tokeep=['TMAX','TMIN','TAVG']
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
    
#All analysis fucntions for a station are done in this object.
#The input is stationdata, which shoudl be self.station_data_claned from the StationLoad object
#the refperiod is a 2 element list of start and end, reference period of analysis
#autorun simply immediately runs the KPI key_metrics function if True 

class StationAnalyzer :
    def __init__(self,stationdata,refperiod=['2020-01-31','2020-12-31'],autorun=True):
        
        #This YAML file contains a great deal of static information, 
        #such as directory information. 
        yaml_file = open("load_stats_static.yaml")
        self.yaml = yaml.load(yaml_file, Loader=yaml.FullLoader)
    
        #Creating the assigned value containing every single element. 
        self.all_data_long=stationdata
        #Then assining off TMID, TMIN, and TMAX. 
        self.tmid_array = self.all_data_long[np.where(self.all_data_long[:,5]==self.yaml['TMID_INDEX'])][:,:5]
        self.tmin_array = self.all_data_long[np.where(self.all_data_long[:,5]==self.yaml['TMIN_INDEX'])][:,:5]
        self.tmax_array = self.all_data_long[np.where(self.all_data_long[:,5]==self.yaml['TMAX_INDEX'])][:,:5]
        
        #This then creates the yearly averages.
        #First, find the beginning and end years, and the number of years total.
        self.alltime_startyear=np.nanmin(self.tmax_array[:,0])
        self.alltime_endyear=np.nanmax(self.tmax_array[:,0])
        self.alltime_number_years=self.alltime_endyear-self.alltime_startyear
        
        #Converting strings to pandas date time
        self.refstart=pd.DatetimeIndex([refperiod[0]])
        self.refend=pd.DatetimeIndex([refperiod[1]])
        
        
        
        
        #Runs the key_metrics function if desired
        if autorun==True:
            self.kpi=self.key_metrics()
            #print(self.kpi)
        
        
    #This function finds the Day of year, assuming every month has 31 days
    #taking in a pd.DatetimeIndex object
    def find_doy(self,pddate):
        doy = pddate.day+(pddate.month-1)*31
        return doy
    
    
    #THIS ISN"T WORKING YET. SOME KIDN OF BUG WHEN TAKING IN THE REFSTART AND REFEND.
    #This function creates beginning and end dates for two different pd datetimeIndex
    #objects, and creates a range of 
    def find_range(self,pddate1,pddate2):
        #First, create simple forms of the Day of year and Year. 
        doy1 = self.find_doy(pddate1)
        doy2 = self.find_doy(pddate2)
        year1=pddate1.year[0]
        year2=pddate2.year[0]
        #In the simple case, the referenc eperiod is in the same year,
        #The array is very simple, jsut each row corresonding to the same year,
        #with a year and DOY.
        if year1==year2:
            if doy1>doy2:
                
                print("Your reference period is improperly defined.")
                print("THe first date must be before the second date.")
                sys.exit("Break Error due to bad dates. .")
            alldays = np.arange(doy1,doy2+1)
            allyears = year1[0]
            returner = np.transpose(np.array([np.repeat(year1,len(alldays)),alldays]))
            
        #If year1 is before year2, the array becomes more compelx. 
        if year1 < year2:
            #First, define the days for the first year. This includes the days from the
            #begining of reference period to the end of the  eyra of teh reference period. 
            year1days=np.arange(doy1,372+1)
            returner = np.transpose(np.array([np.repeat(year1,len(year1days)),year1days]))
            
            #some specical steps have to eb taken, if there are more than 2 adjacent year included.
            if ( year2 - year1 ) > 1 :
                #First, define all the days in the list. 
                dayskarr = np.arange(1,372+1)
                for k in np.arange(year1+1,year2):
                    #Then, over all years in between year 1 and eyar 2, create a new array
                    #This one has DOY frmo 1 to 372 and every yera in the intermediate period. 
                    
                    yearkarr = np.transpose(np.array([np.repeat(k,len(dayskarr)),dayskarr]))
                    returner = np.append(returner,yearkarr,axis=0)
            
            #Then, define the days for the seconds year. This includes the days from the
            #end of reference period to the beginning of the year of the last year of reference. 
            year2days=np.arange(1,doy2+1)
            year2arr= np.transpose(np.array([np.repeat(year2,len(year2days)),year2days]))
            returner = np.append(returner,year2arr,axis=0)
            
        if year1>year2:
            
            print("Your reference period is improperly defined.")
            print("THe first year must be before the second year.")
            sys.exit("Break Error due to bad dates. .")
            
            
                
        return returner    
    
    #This creates a table (pd dataframe) of the key metrics of interest
    def key_metrics(self):       
        #This first does the calculations, finding TMID avreage, TMAx, and TMIn,
        #of the entire record.
        whole_tmid_mean=np.nanmean(self.tmid_array[:,4])
        whole_tmax=np.nanmax(self.tmax_array[:,4])
        whole_tmin=np.nanmin(self.tmin_array[:,4])
        
        #Calcualtes the frequency of extreme temperatures each year. 
        hot_days = np.count_nonzero(self.tmax_array[:,4]>=self.yaml['HOTDAYS'])/(365)
        cold_days = np.count_nonzero(self.tmin_array[:,4]<=self.yaml['COLDDAYS'])/(365)
        
        #This finds the index location of the max temperature
        #First, by finding its index locations. 
        maxloc=np.where(self.tmax_array[:,4]==whole_tmax)
        #Then, creating a list (since there may be more than one) of days of maxium temp. 
        maxyears=   self.tmax_array[maxloc][:,0]
        maxmonths=   self.tmax_array[maxloc][:,1]
        maxdays=   self.tmax_array[maxloc][:,2]
        
        #THen creating a list of strings including every day that was at maximum. 
        maxdate=[]
        for i in np.arange(0,len(maxyears)):
            maxdatestr=str(str(int(maxyears[i]))+'-'+str(int(maxmonths[i]))+"-"+str(int(maxdays[i])))
            maxdate.append(maxdatestr)
            
        #This finds the index location of the min temperature
        #First, by finding its index locations. 
        minloc=np.where(self.tmin_array[:,4]==whole_tmin)
        #Then, creating a list (since there may be more than one) of days of maxium temp. 
        minyears=   self.tmin_array[minloc][:,0]
        minmonths=   self.tmin_array[minloc][:,1]
        mindays=   self.tmin_array[minloc][:,2]
        
        #THen creating a list of strings including every day that was at maximum. 
        mindate=[]
        for i in np.arange(0,len(minyears)):
            mindatestr=str(str(int(minyears[i]))+'-'+str(int(minmonths[i]))+"-"+str(int(mindays[i])))
            mindate.append(mindatestr)
     
        #This creates a dataframe of data points.
        returner = pd.DataFrame(
            {
             "first_year_in_record":[self.alltime_startyear],
             "last_year_in_record":[self.alltime_endyear],
                "TMID_av_F_alltime":[whole_tmid_mean],
              "TMAX_F_alltime":[whole_tmax],
             "TMAX_alltime_date":[maxdate],
             "Annual_days_over_"+str(self.yaml['HOTDAYS']):[hot_days],
             "TMIN_F_alltime":[whole_tmin],
             "TMIN_alltime_date":[mindate],
             "Annual_days_under_"+str(self.yaml['COLDDAYS']):[cold_days],
             
             }
            )
        
        return returner   
    #end key_metrics
     #This creates several charts of interest.
    def key_charts(self):
  
        #Creates a histogram of daily TMID values. 
        plt.hist(self.tmid_array[:,4],bins=30,density=True)
        plt.title('Freqency Histogram of TMID')
        plt.xlabel('Fraction of Total at this TMID')
        plt.ylabel('Temperature, F (TMID)')
        plt.show()
        
        #This creates an array with mean of all months. 
        yearsout=np.empty([0,2])
        for year in np.arange(self.alltime_startyear,self.alltime_endyear+1):
                 
            #Select only the index values wehre the year is equal to year.
            index=np.where(self.tmid_array[:,0]==year)
            #Then, selects data for this year
            #All columns are pulled, since we have to filter to ensure ther eare 
            #12 months for every year evaluated. 
            all_years_data=self.tmid_array[index]
            
            #The data is only written into the new numpy array, if
            #the number of months for the year is 12.
            #This is evaluate dby counting the number of unique days
            #in each year, which must exceed 330. 
            if len(all_years_data)>11*30:
                #And evaluates their mean. 
                all_years_means=np.nanmean(all_years_data[:,4])
                row=np.array([year,all_years_means])
                yearsout=np.append(yearsout,[row],axis=0)
        yearsmean=yearsout
        
        
        #Creates a plot over time of the annual average values. 
        plt.plot(yearsmean[:,0],yearsmean[:,1])
        plt.title('Annual Average of Monthly TMID')
        plt.xlabel('Year')
        plt.ylabel('Annual Average Temperature, F (TMID)')
        plt.show()
        
        #This creates an array with mean of all months. 
        monthout=np.empty([12,2])
        for month in np.arange(1,13):
            
            #Select only the index values wehre the month is equal to month. 
            index=np.where(self.tmid_array[:,1]==month)
            #Then, selects values for this month
            all_months_data=self.tmid_array[index,4]
            #And evaluates their mean. 
            all_months_means=np.nanmean(all_months_data)
            monthout[month-1,0]=month
            monthout[month-1,1]=all_months_means
            
        #Creates a plot over one year of the monthyl average TMID temperatures. 
        plt.plot(monthout[:,0],monthout[:,1])
        plt.title('Monthly Average of TMID')
        plt.xlabel('Month')
        plt.ylabel('Monthly Average Temperature, F (TMID)')
        plt.show()
         
                