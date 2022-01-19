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
#Drop out years with missing data.
returner1 = bzdata.station_data
#First, list out the unique years.
uniqueyears = np.unique(returner1[:,0])
#This defines the minimum number of days of data which msut be present in eveyr month
#for the year to be included.
min_days=15

showme = True

#This initializes the "to keep" list -- the list of data to be included.
tokeep = np.empty([0,1],dtype=int)
#Loop oer all years except the last one.
for year in uniqueyears:
    #As a default, the year is not returned.
    return_this_year = False    
    #Select only the index values wehre the year is equal to year.
    index=np.where(returner1[:,0]==year)
    #Then, selects data for this year
    #All columns are pulled, since we have to filter to ensure ther eare 
    #12 months for every year evaluated. 
    yearsubset=returner1[index]
    
    #Finding all unique months in this year.
    uniquemonths = np.unique(yearsubset[:,1])
    #Coutns the number of months.
    norows = len(uniquemonths)
    #Print an update, if the years is removed.
    if (showme) and (norows < 12) and (year != uniqueyears[-1]):
        print("I dropped "+str(year)+" for having only "+str(norows)+" months of data.")
    #Checks that 12 months are actualyl present in the data.
    if norows == 12:
        #This becomes True for now.
        return_this_year=True
        #Then, checks there are adequate data in each month for this year.
        for month in uniquemonths:
            #Creates a list jsut including tehse months. 
            monthssubset = yearsubset[np.where(yearsubset[:,1]==month)]
            #Finds all numbered values.
            listoftmax = ~np.isnan(monthssubset[:,4])
            listoftmin = ~np.isnan(monthssubset[:,5])
            listoftmid = ~np.isnan(monthssubset[:,6])
            #Counts the number of numbers. 
            tmax_number = np.count_nonzero(listoftmax)
            tmin_number = np.count_nonzero(listoftmax)
            tmid_number = np.count_nonzero(listoftmax)
            #Then checks if each is over the require dminium.
            #Every data point must be present at the right threshold.
            tmax_enough = (tmax_number >= min_days)
            tmin_enough = (tmax_number >= min_days)
            tmid_enough = (tmax_number >= min_days)
            
            if (tmax_enough == False) or (tmin_enough==False) or (tmid_enough ==False):
                return_this_year = False
                if showme: print("I dropped "+str(year)+" for having less than "+str(min_days)+" days in month #"+str(month))
    if (return_this_year ==True) or (year == uniqueyears[-1]):
        #If everything is true; then add it to the list of rows to keep. 
        tokeep=np.append(tokeep,index)
returner1 =returner1[tokeep]  
newuniqueyears=np.unique(returner1[:,0])