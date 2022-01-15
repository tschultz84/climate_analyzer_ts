# -*- coding: utf-8 -*-
"""
Created on Wed Dec 29 13:10:53 2021

@author: 14154
"""
import Station_Loader_ts as load
import Station_analyzer_ts as analyze


#%%
bzdata=load.LoadStation([45.647256643331126,-111.04060494981753],10,True,True) #Loading in Bozeman, MT coordinates

#%%

date1 = '2021-6-1'
date2 = '2021-8-30'
bzcalc=analyze.StationAnalyzer(bzdata.station_data_clean,date1,date2,display=True)

#test=bzcalc.tmid_selection(bzcalc.refperiod)
#print(test)


#%%
obdata=load.LoadStation([37.755663644,-122.506497974],10,True,True) #Loading in Ocean Beach, SF coordinates
#obdata.calculate_tmid_new(obdata.all_data_np)


#%%
date1 = '2010-1-1'
date2 = '2020-12-1'
obcalc=analyze.StationAnalyzer(obdata.station_data_clean,date1,date2,display=True)

#print(obcalc.kpi)
#obcalc.key_charts()

#%%
sddata=load.LoadStation([32.741947,-117.239571],10,True,True) #Loading in Ocean Beach, San Diego coordinates
#%%
date1 = '2000-1-1'
date2 = '2020-12-31'
sdcalc=analyze.StationAnalyzer(sddata.station_data_clean)

#%%

dvtdata=load.LoadStation([32.659167, -116.099167],10,True,True) #Loading in Desert View Tower data
#%%
date1 = '2000-1-1'
date2 = '2020-12-31'
dvtcalc=analyze.StationAnalyzer(dvtdata.station_data_clean,date1,date2,display=True)
#print(dvtcalc.kpi)
#dvtcalc.key_charts()
#print(dvtcalc.maxdate)

#

#%%



