# -*- coding: utf-8 -*-
"""
Created on Wed Dec 29 13:10:53 2021

@author: 14154
"""



#%%
bzdata=LoadStation([45.647256643331126,-111.04060494981753],10,True,True) #Loading in Bozeman, MT coordinates

#%%

date1 = '2020-1-1'
date2 = '2020-12-31'
bzcalc=StationAnalyzer(bzdata.station_data_clean,date1,date2)

#test=bzcalc.tmid_selection(bzcalc.refperiod)
#print(test)


#%%
obdata=LoadStation([37.755663644,-122.506497974],10,True,True) #Loading in Ocean Beach, SF coordinates
#obdata.calculate_tmid_new(obdata.all_data_np)


#%%
date1 = '2010-1-1'
date2 = '2020-12-31'
obcalc=StationAnalyzer(obdata.station_data_clean,date1,date2)

#print(obcalc.kpi)
#obcalc.key_charts()

#%%
sddata=LoadStation([32.741947,-117.239571],10,True,True) #Loading in Ocean Beach, San Diego coordinates
#%%
date1 = '2015-1-1'
date2 = '2020-12-31'
sdcalc=StationAnalyzer(sddata.station_data_clean)

#%%

dvtdata=LoadStation([32.659167, -116.099167],10,True,True) #Loading in Desert View Tower data
#%%
dvtcalc=StationAnalyzer(dvtdata.station_data_clean)
#print(dvtcalc.kpi)
#dvtcalc.key_charts()
#print(dvtcalc.maxdate)

#

#%%



