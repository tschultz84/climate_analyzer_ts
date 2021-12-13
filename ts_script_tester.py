# -*- coding: utf-8 -*-
"""
Created on Mon Dec 13 13:07:09 2021

@author: 14154
"""


import shapely
from shapely.geometry import Point,MultiPoint
from shapely.ops import nearest_points
import geopandas as gpd
import pandas as pd
#%%
#DATAPATH is just where ghcnd-stations and ghcnd-inventory are located on your hard drive.
stDATAPATH = "C:\\ts_big_data_files\\"
test=pd.read_csv(stDATAPATH+'ghcnd_station_master_ts.csv')

#%%

testlatlon =[ 45.676998,-111.042931] #Lat Lon of Bozeman Mt