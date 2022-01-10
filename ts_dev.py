# -*- coding: utf-8 -*-
"""
Created on Thu Dec 23 06:51:09 2021

@author: 14154
"""



#%%

bztmiddata=bzcalc.tmid_np


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
data1=bzcalc.tmid_ref_data[:,4]
data2=bzcalc.tmid_base_data[:,4]

data1=data1[~np.isnan(data1)]
data2=data2[~np.isnan(data2)]

print(np.var(data2)/np.var(data1))

print(stats.ttest_ind(a=data1, b=data2, equal_var=True))

    
