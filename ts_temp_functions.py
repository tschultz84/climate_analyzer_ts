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


#%%
#This function takes as an input, a Point(lon,lat) that is the 
#reference point of interest
#It then searches for the nearest point in our list of points,
#searching distance by miles distance
#It returns the top 10 nearest stations, ordered
#This is returned in the foramt of a dataframe.
#The first row of the returned dataframe is the reference point. 

def nearest_station(pt):

   #DATAPATH is just where ghcnd-stations and ghcnd-inventory are located on your hard drive.
   stDATAPATH = "C:\\ts_big_data_files\\"
   #This initializes the actual data to search, from the list of stations
   #It only reads in the columns which are necessary to search by distance. . 
   df = pd.read_csv(stDATAPATH+'ghcnd_station_master_ts.csv')[['ID','Name','Latitude','Longitude']].drop_duplicates()
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
   

   return refdf.append(returner)

#%%
#This function returns the miles between two poitns
#def ts_latlon_distance(lat1, lat2, lon1, lon2):
def ts_latlon_distance(latlon1,latlon2):
    # radians which converts from degrees to radians.
    lon1 = radians(latlon1[1])
    lon2 = radians(latlon2[1])
    lat1 = radians(latlon1[0])
    lat2 = radians(latlon2[0])
    
    """
    # The math module contains a function named
    # radians which converts from degrees to radians.
    lon1 = radians(lon1)
    lon2 = radians(lon2)
    lat1 = radians(lat1)
    lat2 = radians(lat2)
      """
    # Haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat / 2)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2)**2
 
    c = 2 * asin(sqrt(a))
    
    # Radius of earth in kilometers. Use 3956 for miles
    r = 3956
      
    # calculate the result, returned in miles
    return(c * r)
     
