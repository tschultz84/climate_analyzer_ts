# -*- coding: utf-8 -*-
"""
Created on Thu Dec 23 06:51:09 2021

@author: 14154
"""
#%%
test10=LoadStation([45.647256643331126,-111.04060494981753],10,True,True) #Loading in Bozeman, MT coordinates
print(test10.station_data)
print(test10.station_data_clean)
#%%

data=test10.station_data_clean


#%%
data11=data1[data1['Year']==1952]
print(data11[data11['Element']=='TMAX'])
print(data11[data11['Element']=='TMIN'])
print(data11[data11['Element']=='TMID'])
