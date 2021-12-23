# -*- coding: utf-8 -*-
"""
Created on Thu Dec 23 06:51:09 2021

@author: 14154
"""
import time
#%%
#This function takes as an input, a Point(lon,lat) that is the 
#reference point of interest
#It then searches for the nearest point in our list of points,
#searching distance by miles distance
#It returns the top 10 nearest stations, ordered
#This is returned in the foramt of a dataframe.
#The first row of the returned dataframe is the reference point. 

def nearest_station_dev(pt):
   #startTime = time.time()
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
   
   print("Searching closest station among "+str(len(df))+" stations.")
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
       d1er.append(ts_latlon_distance([pt.y,pt.x],otw.iloc[i][['Latitude','Longitude']]))
      
   #This appends the data series to the final return dataframe
   otw['Miles_from_Ref']=d1er
   #and then sorts by miles from ref, to find the closest stations
   returner=otw.sort_values(by='Miles_from_Ref')
   
   #Then return the closest 10 stations.
   returner=returner[0:10]
   
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
   return refdf.append(returner)


#%%
test = Point(-111.04060494981753,45.647256643331126)
rf=nearest_station_dev(test)
#nearest_station_dev_1(test)