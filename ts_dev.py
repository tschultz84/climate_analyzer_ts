# -*- coding: utf-8 -*-
"""
Created on Thu Dec 23 06:51:09 2021

@author: 14154
"""
import os
import numpy as np

#%%
mttest = Point(-111.04060494981753,45.647256643331126) #Bozeman lat lon
sdtest = Point(-117.239571,32.741947) #Ocean Beach SD lat lon
sftest = Point(-122.506497974,37.755663644  )
#rf=nearest_station(test)
#nearest_station_dev_1(test)

#%%
test1 = FindStation(sftest,True)
#%%
import os
import sys
import yaml
#%%
#This object contains the functions to read in key data about a station.
#Its input is a station ID. 
class StationReader:
    def __init__(self,stationID):
        self.stid = stationID
        
        #This YAML file contains a great deal of static information, 
        #such as directory information. 
        yaml_file = open("load_stats_static.yaml")
        self.yaml = yaml.load(yaml_file, Loader=yaml.FullLoader)
        
        
        self.datapath= self.yaml['DATAPATH']

        
    #This loads the station, using the station ID which is passed.
    def load_station(self,station):
        
        startTime = time.time()
        if station == None:
            station=self.stid
        #This assigns the filename where station info is located.
        filename=self.datapath+str(station)+'.dly'
        #This checks the file exists and breask if it does not.
        #(If it doesn't exist, the submitted station ID was in error.)
        if os.path.exists(filename) == False:
            print("Bad Station ID. The file called "+filename+" does not exist.")
            print("Please check your station ID and re-submit.")
            sys.exit("Break Error in load_station of StationReader: Bad Station ID.")
        #You have designed this script so it only reads in the lines where data is useful
        #and excludes other elements. 
        #When testing reading in the very largest data file, you would typically save
        #0.02 to 0.06 seconds during loading. Worth every bit...
        
        #This defines the pre-loading data, i.e. which columns are used to skip rows.
        shortcols = [ (17, 21)]
        shortnames =  np.array([ 'Element'])
        
        #This is to test whether the prelim, read-in step, helps or hinder the read time
        
        print("Stripping out extraneous rows before reading everything in...")
        #This then reads the file, but only the columns which generate the rows to skip. 
        prelim = pd.DataFrame()    
        prelim=pd.read_fwf(filename,colspecs=shortcols)
        prelim.columns=shortnames
        #This generates the rows to skip.   
        #Currently you are keeping only TMAX and TMIN
        tokeep=['TMAX','TMIN']
        skiprows1 = prelim[(~prelim['Element'].isin(tokeep))].index
        
        
        #Then, the data is read in.
        inter0 = pd.DataFrame()      
        #This then Reads out the next file.  
        inter0 = pd.read_fwf(filename,colspecs=self.yaml['datacolnums'],skiprows=skiprows1+1)
        inter0.columns=self.yaml['datacolnames']
            
 
            
        self.station_data=inter0
        executionTime = (time.time() - startTime)
        print('Time taken to load station data using load_station: ' + str(executionTime))
        return inter0

       
#%%
test5=StationReader('AU000005901')
test6=test5.load_station(None)
#%%
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
