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
bzdata=LoadStation([45.647256643331126,-111.04060494981753],10,True,True) #Loading in Bozeman, MT coordinates

#%%
bzcalc=StationAnalyzer(bzdata.station_data_clean)
print(bzcalc.kpi)
bzcalc.key_charts()

#%%
obdata=LoadStation([37.755663644,-122.506497974],10,True,True) #Loading in Ocean Beach, SF coordinates
#%%
obcalc=StationAnalyzer(obdata.station_data_clean)

print(obcalc.kpi)
obcalc.key_charts()

#%%
test12=LoadStation([32.741947,-117.239571],10,True,True) #Loading in Ocean Beach, San Diego coordinates
#%%
dvtdata=LoadStation([32.659167, -116.099167],10,True,True) #Loading in Desert View Tower data
#%%
dvtcalc=StationAnalyzer(dvtdata.station_data_clean)
print(dvtcalc.kpi)
#testdvt.key_charts()

#%%

#%%



