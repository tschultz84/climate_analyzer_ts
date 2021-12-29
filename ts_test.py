# -*- coding: utf-8 -*-
"""
Created on Wed Dec 29 13:10:53 2021

@author: 14154
"""

#%%
from shapely.geometry import Point,MultiPoint,Polygon
mttest = Point(-111.04060494981753,45.647256643331126) #Bozeman lat lon
sdtest = Point(-117.239571,32.741947) #Ocean Beach SD lat lon
sftest = Point(-122.506497974,37.755663644  )
#rf=nearest_station(test)
#nearest_station_dev_1(test)

#%%
test1 = FindStation(mttest,False)
test10=test1.nearest_station()

test55=StationReader(test1.id_closest_st)
test66=test55.load_station(None)
#%%
tes7=test5.load_station('CA1AB000158')

#%%
test10=LoadStation([45.647256643331126,-111.04060494981753],10,True,True) #Loading in Bozeman, MT coordinates
print(test10.closest_stations)
print(test10.station_data)

#%%
test11=LoadStation([37.755663644,-122.506497974],20,True,True) #Loading in Ocean Beach, SF coordinates

print(test11.station_data)
