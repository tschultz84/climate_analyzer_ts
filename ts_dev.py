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
bzpand = bzdata.all_data_1

allyears=np.unique(bzpand[:,0])
alleles=np.unique(bzpand[:,2])

for k in alleles:
    for i in allyears:
        print(k)
        print(i)
        
#%%
tmax=bzdata.tmaxval
tmin=bzdata.tminval

tmax[0,0:4]==tmin[0,0:4]   
outer= []
for i in np.arange(0,len(tmax)):
    year=tmax[i,0]
    doy = tmax[i,3]
    index = np.where((tmin[:,0]==year)&(tmin[:,3]==doy))
    if np.shape(index)[1]>0:
        tminval = tmin[index,4][0][0]
        outer.append([tmax[i,0],tmax[i,1],tmax[i,2],tmax[i,3],tmax[i,4],tminval])
outer=np.asarray(outer)
means = np.nanmean(np.array([outer[:,4],outer[:,5]]),axis=0)
dates=outer[:,0:4]
returner1 = np.insert(dates,4,means,axis=1)

#%%
start=bzcalc.refstart
end=bzcalc.refend

start=pd.DatetimeIndex(['2018-12-31'])

#This function finds the Day of year, assuming every month has 31 days
#taking in a pd.DatetimeIndex object
def find_doy(pddate):
    doy = pddate.day+(pddate.month-1)*31
    return doy[0]

print(find_doy(start))

#%%

date1=pd.DatetimeIndex(['2020-12-30'])
date2=pd.DatetimeIndex(['2019-1-11'])


#This function creates beginning and end dates for two different pd datetimeIndex
#objects, and creates a range of 
def find_range(pddate1,pddate2):
    #First, create simple forms of the Day of year and Year. 
    doy1 = find_doy(pddate1)
    doy2 = find_doy(pddate2)
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
        
        
print(find_range(date1,date2)[350:400])

    
