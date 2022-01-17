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
# String with URL:
url_dir = 'https://www.ncei.noaa.gov/data/global-historical-climatology-network-daily/access/'
url_name = 'ACW00011604.csv'
columns=['DATE',"TMAX","TMIN","TAVG"]
# First example to read csv from URL
#df = pd.read_csv(yaml['NOAA_URL']+url_name,usecols=columns)
df = pd.read_csv(yaml['NOAA_URL']+url_name,usecols=columns)
#df.columns()
#%%
import requests
response = requests.get(url_dir+'XXX')
#%%
yaml_dir="C:\\Users\\14154\\OneDrive\\Python\\climate_mapper\\python\\climate_analyzer_ts\\"
yaml_file = open(yaml_dir+"load_stats_static.yaml")
yaml = yaml.load(yaml_file, Loader=yaml.FullLoader)