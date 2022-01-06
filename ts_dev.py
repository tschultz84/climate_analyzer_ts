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


#%%

bztmiddata=bzcalc.tmid_np
#%%
inter=[]
for k in np.arange(0,len(bztmiddata)):
    
    year=int(bztmiddata[k,0])
    month=int(bztmiddata[k,1])
    monthdata=np.transpose(bztmiddata[k,2:])
    
    for i in np.arange(1,31):
        doy = i + (month-1)*31
        inter.append([year,month,i,doy,monthdata[i-1]])
        
    
inter=np.asarray(inter)

#%%
ele='TMAX'
if(ele == "TMIN"):
    eler=0
if(ele == "TMAX"):
    eler=1
if(ele == "TMID"):
    eler=2