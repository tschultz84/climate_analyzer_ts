# -*- coding: utf-8 -*-

import pandas as pd
import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
from sklearn.linear_model import LinearRegression
from scipy.stats.stats import pearsonr
from IPython.display import display_html 


#All analysis fucntions for a station are done in this object.
#The input is stationobj, which is an object output from the Station_Loader_ts object
#the refst and refend are both strings, marking the beginning and end of the reference period of analysis
#alpha is the alpha value (p-value threshold determining statistical signifiance)

#display prints the function outputs
import json
import os
import Station_Loader_ts as load
import imp
imp.reload(load)

class StationAnalyzer :
    def __init__(self,point,display=True,alpha=0.01):
        #Load the weather station. 
        stationobj=load.LoadStation(point)

        #Initialize some variables.
        self.alpha=alpha
        self.display=display
        self.stationobj = stationobj
        
        #Creating the assigned value containing every single element. 
        self.all_data_long=stationobj.station_data
        #Then assigning off TMID, TMIN, and TMAX values, which are different
        #columns of the input data.
        self.tmid_array = self.all_data_long[:,[0,1,2,3,6,7]]
        self.tmin_array = self.all_data_long[:,[0,1,2,3,5,7]]
        self.tmax_array = self.all_data_long[:,[0,1,2,3,4,7]]
        
        #Initializes the reference dates. It initializes the reference period to 2015-1-1 to 2020-12-31.
        self.change_ref_dates([2015,2020])
        
        #Runs the key_metrics function if desired.
        #if self.display==True:
            #Generate variables. 
            #self.key_metrics()
            #show everything in a nice format. 
            #display_html(self.key_metrics_table)
            #display_html(self.key_stats)
            #self.key_charts()
                
    #BEGIN FUNCTIONS
    """
    DATA CLEANING AND WRANGLING FUNCTIONS
    The following set of functions make various transformations on data
    and create new arrays in order to complete various analyses.
    """  
    #Change Reference Period Dates.
    #It takes in year_series, a list of [year_start,year_end] which is the beginning and end of the period.
    #then date_series, list of [[month, day],[month, day]] which is the beginning and end of a subannual reference.
    #If date_series != list, then it assumes the beginning and end are the whole year.
    def change_ref_dates(self,year_series,date_series=-1):
        #If date_series is not a list, then set the date differently.
        if type(date_series)!=list:
            date_series = [[1,1],[12,31]] #Set to the entire year.
        
        #Create the dates in date format.
        self.ref_date_1=pd.to_datetime(f"{year_series[0]}-{date_series[0][0]}-{date_series[0][1]}")
        self.ref_date_2=pd.to_datetime(f"{year_series[1]}-{date_series[1][0]}-{date_series[1][1]}")

        #This simply throws an error if the dates are not ordered properly.
        if self.ref_date_1.to_julian_date()>self.ref_date_2.to_julian_date():        
            print("Your period dates are improperly defined. The first date must be before the second date.")
            sys.exit("Break Error due to bad dates.")
        print(f"Reference dates have been re-defined.")
        print(f"Years: {self.ref_date_1.year} to {self.ref_date_2.year}.")
        print(f"Period Subset Included: {self.ref_date_1.month}-{self.ref_date_1.day} to {self.ref_date_2.month}-{self.ref_date_2.day}.")
        
        #Regenerate the new periods of interest. 
        self.refperiod=self.find_time_range(ref=True) #Creating the arrayed reference period. 
        self.baseperiod = self.find_time_range(ref=False)  #creates baseline range.

        #Create a series summarizing information.
        self.period_information=pd.Series(data={
            "Reference Period Start":self.ref_date_1,
            "Reference Period End":self.ref_date_2,
            "Reference Period Start (Julian)":self.ref_date_1.to_julian_date(),
            "Reference Period End (Julian)":self.ref_date_2.to_julian_date(),
            "First Unique DOY in Reference Period":self.ref_date_1.day_of_year,
            "Last Unique DOY in Reference Period":self.ref_date_2.day_of_year
        })
        
        #Creates the actual values of TMID in the periods.
        self.tmid_ref_data=self.tmix_selection(self.refperiod,self.tmid_array)
        self.tmid_base_data=self.tmix_selection(self.baseperiod,self.tmid_array)
        self.all_years_mean=self.alltime_ref_annual_averages()
        self.all_months_mean = self.alltime_annual_monthly_averages()
        print(f"I re-generated the baseline and reference period datasets accordingly. Their respective shapes: {self.tmid_ref_data.shape} and {self.tmid_base_data.shape}.")

    #This function creates beginning and end dates for two different pd datetimeIndex
    #objects, and creates a range of days, in the form of [Julian_date_1,Julian_date_2, etc...]
    #It automatically does this based on the the baseline and reference periods. 
    #The 2nd argument, reference, just defines if the range is in the baseline or not.
    def find_time_range(self,ref=True):

        if ref==True: #If the reference period, use the reference period to define dates.
            year_1,year_2 = self.ref_date_1.year,self.ref_date_2.year
            
        if ref!=True: #If not True, then use the baseline period.
            begdate= self.stationobj.station_information['Earliest year in station record']
            enddate = self.stationobj.station_filters["Last Year of Baseline"]
            year_1,year_2=pd.to_datetime(begdate,format="%Y").year,pd.to_datetime(enddate,format="%Y").year
        #The days of year are the same in both analyses
        doy_1, doy_2 = self.ref_date_1.day_of_year, self.ref_date_2.day_of_year
            
        #Create the output array by looping over the whole reference period.
        return_arr = np.empty([0,1]) #Create an initialize array. 
        for year in np.arange(year_1,year_2+1): #Loop over all years.
            for day in np.arange(doy_1,doy_2+1): #Loop over all Days of Year. 
                arr=np.array([pd.to_datetime(f"{int(year)}-{int(day)}",format="%Y-%j").to_julian_date()])
                return_arr = np.append(return_arr,[arr],axis=0)
        
        return return_arr 
        
    #This function selects the days corresponding to "range1" from the dataset "data".
    #range1 is an array in the format of self.refperiod, a long list of rows in the form [julian_date]
    #data is  data, in the format of self.tmid_array
    def tmix_selection(self,range1,data):
        #Rename the input variable, which are subject to change.
        out=data
        #Now use Julian dates to filter down to the correct range.
        index = np.in1d(out[:,5],range1)
        #Select only the subset of time defined by index.
        out=out[index]
        return out
            
    """ANALYSIS FUNCTIONS"""
    """
    all_time_ref_annual_averages creates an array of average temperatures across the entire time period.
    The average is evaluated only over the days provided in the reference period.
    So if you provide only a month's worth of data, you see only the averages for that
    month.'
    """
    def alltime_ref_annual_averages(self,):
        yearsout=np.empty([0,2])
        #Define the list of DOY to be included.
        doy_list=np.arange(
            self.period_information["First Unique DOY in Reference Period"],
            self.period_information["Last Unique DOY in Reference Period"]+1)
        #This extracts from tmid_array only the doy present in the reference period.
        indexer = np.in1d(self.tmid_array[:,3],doy_list)
        yearsdata=self.tmid_array[indexer]
        
        uniqueyears=np.unique(yearsdata[:,0])
        for year in uniqueyears[:-1]:   
            #Select only the index values wehre the year is equal to year.
            index=np.where(yearsdata[:,0]==year) 
            all_years_data=yearsdata[index]
            
            #And then evaluates their mean. 
            all_years_means=np.nanmean(all_years_data[:,4])
            row=np.array([year,all_years_means])
            yearsout=np.append(yearsout,[row],axis=0)
        return yearsout
    
    #This creates an array with monthyl average temperatures. 
    def alltime_annual_monthly_averages(self):
        monthout=np.empty([12,2])
        for month in np.arange(1,13):
            #Select only the index values wehre the month is equal to month. 
            index=np.where(self.tmid_array[:,1]==month)
            #Then, selects values for this month
            all_months_data=self.tmid_array[index,4]
            #And evaluates their mean. 
            all_months_means=np.nanmean(all_months_data)
            monthout[month-1,0]=month
            monthout[month-1,1]=all_months_means
        return monthout

    """STATISTICAL TESTS ----------
    The following functions perform various statistical tests."""
    #This completes a t-test of the reference vs. baseline period TMID values. 
    #This function reviews the TMID data for the reference and baseline,
    #compares them, does a t-test using stats.ttest_ind.
   
    def do_t_test(self):
        #Define the datasets. 
        data1=self.tmid_ref_data[:,4]
        data2=self.tmid_base_data[:,4]
        #Scrub the not as number values. 
        data1=data1[~np.isnan(data1)]
        data2=data2[~np.isnan(data2)]
        
        """
        if standard=True, perform a standard independent 2 sample t-test that 
        assumes equal population variances. 
        If False, perform Welch’s t-test, which does not
        assume equal population variances. This is True by default.
        Before we perform the test, we need to decide
        if we’ll assume the two populations have equal
        variances or not. As a rule of thumb, we can 
        assume the populations have equal variances if 
        the ratio of the larger sample variance to the smaller sample variance is less than 4:1. 
        """
        standard=True
        if (np.var(data2)/np.var(data1) > 4 ) or (np.var(data2)/np.var(data1) <0.25 ):
            standard=False
        result = stats.ttest_ind(a=data1, b=data2, equal_var=standard)
        return result
    
    #This creates a trend line.
    #It takes as input a numpy array containing yearly averages, in the format of self.all_years_mean
    #Its output is a tuple, with the first two being the slope and intercept.
    #and the third element is an array of the year and value of the line according to this regression.
    #If array=False, then this array is not generated or returned 
    def create_trend_line(self,yearsdata,array=True):
        #If array was set improperly, automatically change it to a Boolean true value. 
        if (array != False) & (array != True):
            array = True
        #First, format, and take only the recent years of data.
        recentyears=30 #Using 30 years for the trend.
        x = yearsdata[-recentyears:,0].reshape(-1,1)
        y = yearsdata[-recentyears:,1]
        #This statement creates the variable model as the instance of LinearRegression. 
        model = LinearRegression()
        # call .fit() on model:
        model.fit(x,y) 
        intercept = model.intercept_
        slope = model.coef_
        if array == True:
            return_arr = np.empty([0,2])
            for k in np.arange(0,len(x)):
                year = x[k]
                value = slope*year + intercept
                arr = [year,value]
                return_arr = np.append(return_arr,arr)
            return_arr = np.reshape(return_arr,[-1,2])
            
            return slope,intercept,return_arr
        if array == False:
            return slope,intercept
    
    """METRIC ANALYSES FUNCTIONS
    These functions return various metrics of analyses used in key_metrics.
    """
    
    #This function takes an input station_data in the format of tmid_array (or tmax_array)
    #and returns the dates when the maximum or minimum temperature occured in the timeframe..
    #If MAX=TRUE, then it finds the date of the temperature maximum in the range of dates in tmix_array.
    #If MAX=FALSE, then it does this for the temperature minimum
    def find_max_min_info(self,tmix_array,MAX=True):
       #First, find tmax or tmin value.
        if MAX == True:
            value = np.nanmax(tmix_array[:,4])
        if MAX == False:
            value = np.nanmin(tmix_array[:,4])
        
        #This finds the index location of the max/min temperature
        #First, by finding its index locations. 
        value_loc=np.where(tmix_array[:,4]==value)
        #Then, creating a list (since there may be more than one) of days of max/min temp. 
        valyears=   tmix_array[value_loc][:,0]
        valmonths=  tmix_array[value_loc][:,1]
        valdays=   tmix_array[value_loc][:,2]
        
        #THen creating a list of strings including every day that was at the value. 
        valuedates=[]
        for i in np.arange(0,len(valyears)):
            valdatestr=str(str(int(valyears[i]))+'-'+str(int(valmonths[i]))+"-"+str(int(valdays[i])))
            valuedates.append(valdatestr)
        return valuedates
    #PICK UP CLEANUP HERE
    #This creates a table (pd dataframe) of the key metrics of interest
    def key_metrics(self):       
        #ALL TIME METRICS -- first finds the metrics over the whole period.
        #This first does the calculations, finding TMID avreage, TMAx, and TMIn,
        #of the entire record.
        whole_tmid_mean=np.nanmean(self.tmid_array[:,4])
        whole_tmax=np.nanmax(self.tmax_array[:,4])
        whole_tmin=np.nanmin(self.tmin_array[:,4])
        whole_tmid_var=np.nanvar(self.tmid_array[:,4])
        
        #Finds the dates of the maximum and minimum temperatures.  
        maxdate = self.find_max_min_info(self.tmax_array)
        mindate = self.find_max_min_info(self.tmin_array,MAX=False)
       
        #REFERENCE AND BASELINE PERIOD METRICS OF INTEREST
        #First, select the data. 
        refmiddata = self.tmix_selection(self.refperiod, self.tmid_array)
        bmiddata = self.tmix_selection(self.baseperiod, self.tmid_array)
        
        #This then calculates the mean over the ref and baseline periods. 
        ref_tmid_mean=np.nanmean(refmiddata[:,4])
        base_tmid_mean=np.nanmean(bmiddata[:,4])
        
        #This creates variables which can then be extracted. 
        self.refmean = ref_tmid_mean
        self.basemean = base_tmid_mean
                 
        #And variance. 
        ref_tmid_var=np.nanvar(refmiddata[:,4])
        base_tmid_var=np.nanvar(bmiddata[:,4])
        
        #And max, min values in the REF and BASE.
        #First, select the data. 
        refmaxdata = self.tmix_selection(self.refperiod, self.tmax_array)
        refmindata = self.tmix_selection(self.refperiod, self.tmin_array)
        
        bmaxdata = self.tmix_selection(self.baseperiod, self.tmax_array)
        bmindata = self.tmix_selection(self.baseperiod, self.tmin_array)
        #Then, calcualte max and min. 
        ref_max=np.nanmax(refmaxdata[:,4])
        base_max=np.nanmax(bmaxdata[:,4])
        ref_min=np.nanmin(refmindata[:,4])
        base_min=np.nanmin(bmindata[:,4])
        #Then, find max and min dates.
        rmaxdate = self.find_max_min_info(refmaxdata)
        rmindate = self.find_max_min_info(refmindata,MAX=False)
        bmaxdate = self.find_max_min_info(bmaxdata)
        bmindate = self.find_max_min_info(bmindata,MAX=False)
        

        #Creates strings describing the range of dates.
        alltimestr1 = f"{self.tmid_array[-1][0]}-{self.tmid_array[-1][1]}-{self.tmid_array[-1][2]}"
        alltimestr2 = f"{self.tmid_array[-1][0]}-{self.tmid_array[-1][1]}-{self.tmid_array[-1][2]}"
        
        #Creates strings describing the range of dates. 
        firstdayr,lastdayr=self.tmid_ref_data[0][0:3],self.tmid_ref_data[-1][0:3] #Find the year and days of year of the first and last entry in Reference period..
        self.refstring = f"{firstdayr[1]}-{firstdayr[2]} to {lastdayr[1]}-{lastdayr[2]} in {firstdayr[0]} to {lastdayr[0]}" #Generate the string. 
        firstdayb,lastdayb=self.tmid_base_data[0][0:3],self.tmid_base_data[-1][0:3] #Find the year and days of year of the first and last entry in baseline period.
        basestring = f"{firstdayb[1]}-{lastdayb[2]} to {lastdayb[1]}-{lastdayb[2]} in {firstdayb[0]} to {lastdayb[0]}"#Generate the string. 

        #Just pull out the days. 
        #self.refdays = refstr1[5:]+" to "+refstr2[5:]
        #This creates a dataframe of data points.
        returner = pd.DataFrame(
            {
             "Period of Measure":
                 ["All Time: " +alltimestr1+" to "+alltimestr2,
                 "Reference: "+self.refstring,
                 "Baseline: "+basestring],

                "Average TMID over period (F)":[whole_tmid_mean,ref_tmid_mean,base_tmid_mean],
                "Variance of TMID over period (F)":[whole_tmid_var,ref_tmid_var,base_tmid_var],
                
                "Maximum Temperature over period (F)":[whole_tmax,ref_max,base_max],
                "Date of Maximum Temperature":[maxdate,rmaxdate,bmaxdate],
                "Minimum Temperature over period (F)":[whole_tmin,ref_min,base_min],
                "Date of Minimum Temperature":[mindate,rmindate,bmindate],      
             }
            )
        
        ##THese are key values to return out of the function.
        #You can grab these via the object handle. 
        self.alltimestr =alltimestr1+" to "+alltimestr2
        self.whole_mean = whole_tmid_mean
        self.whole_max_info = [whole_tmax,maxdate]
        self.whole_min_info = [ whole_tmin,mindate]
        self.ref_max_info = [ref_max,rmaxdate]
        self.ref_min_info = [ref_min,rmindate]
        
        
        t_test = self.do_t_test()
        #print("T Test Finding:")
        #print(t_test)
        pvalue = t_test[1]
        
        #First, create a string pulling out the difference between reference and baseline.
        self.ref_base_delta = str(int((ref_tmid_mean-base_tmid_mean)*100)/100)
        #And the same with the pvalue.
        self.ref_pvalue = pvalue
        
        
       # print("The Ref period is "+str(int((ref_tmid_mean-base_tmid_mean)*100)/100)+"F warmer than your base period.")
       # print("Is this difference statistically different, ")
        #print(" with a probability that this occured by chance less than "+str(self.alpha))
        #print(" i.e., alpha is "+str(self.alpha)+str("?"))
     
        if pvalue >= self.alpha:
           #print("No.")
           self.ref_stat_sig = False
        if pvalue < self.alpha:
           #print("Yes.") 
           self.ref_stat_sig = True

        trend_data = self.create_trend_line(self.all_years_mean,False)
        
        #This creates a variable to pull out the warming, in F per decade.
        self.trend_data_str=str(round(trend_data[0][0]*10,2))+" Farenheit per decade"
        #This calculates the correlation coefficient and accompanying p-value
        #which states whether this correlation is statistically significant. 
        
        #First, format, and take only the recent years of data.
        recentyears=30 #using 30 years in the trend.
        x5 = self.all_years_mean[-recentyears:,0]
        y5 = self.all_years_mean[-recentyears:,1]
        
        pearsond1 = pearsonr(x5,y5)
                
       # print("Outcome of Pearson Correlation Coefficient Analysis: ")
       # print("Coefficient of Correlation (R): "+str(pearsond1[0]))
       # print("Accompanying P value: "+str(100*pearsond1[1]))
        if pearsond1[1] < self.alpha:
            is_real = True
        if pearsond1[1] >= self.alpha:
            is_real = False
       # print("Is this trened real using Pearson correlation analysis? "+str(is_real))
        #Variables to pull out.
        self.trend_p = pearsond1[1]
        self.trend_is_real = is_real
        #This creates a dataframe of statistical information points.
        key_stats = pd.DataFrame(
            {
                "Metric":
                    ["Reference Minus BasEline Temperature Change",
                     str(30)+ " year warming trend"
                     ],
             "Value":
                 [str(int((ref_tmid_mean-base_tmid_mean)*100)/100)+ "F",
                 self.trend_data_str ],
             "P-Value of Difference":
                 [self.ref_pvalue,
                  self.trend_p],
                 "Statistically Significant?":
                     [self.ref_stat_sig,
                      self.trend_is_real],
                     "Alpha Value":
                         [self.alpha,
                          self.alpha]
                
                }
        )
        #print(key_stats)
        
        self.key_metrics_table =returner
        self.key_stats = key_stats
        return returner   
    
    #end key_metrics
              
     #This creates several charts of interest.
    def key_charts(self):
        plt.rcParams['figure.figsize'] = [8, 6]
        
        #Creates a histogram of daily TMID values. 
        ref_hist=self.tmid_ref_data[:,4]
        base_hist=self.tmid_base_data[:,4]
        plt.hist(ref_hist,bins=30,density=True,label="Reference Period", alpha=0.5)
        plt.hist(base_hist,bins=30,density=True,label="Base Period", alpha=0.5)
        plt.title('Freqency Histogram of Daily TMID over Reference and Base Period')
        plt.xlabel('Fraction of Total at this TMID')
        plt.ylabel('Temperature, F (TMID)')
        plt.legend()
        plt.show()
        
            
        #Creates a plot over one year of the monthyl average TMID temperatures. 
        plt.plot(self.all_months_mean[:,0],self.all_months_mean[:,1])
        plt.title('Monthly Average of TMID- All Time Record')
        plt.xlabel('Month')
        plt.ylabel('Monthly Average Temperature, F (TMID)')
        plt.show()
         
        #This creates a lienar trend to superimpose on the chart.
        trenddat=self.create_trend_line(self.all_years_mean)
        #Pulls out the slope for chart labeling. 
        slope=trenddat[0][0]
        #Then creates the trend line data. 
        trendline=trenddat[2]
        trend_x = trendline[:,0]
        trend_y = trendline[:,1]
        
        #Creates a plot over time of the  average values. 
        #This is averaged over all_years_mean, which only includes days in the 
        #reference period, although for all years.
        nodays1 = len(np.unique(self.refperiod[:,1]))
        if nodays1 > 364 :
            text = "Entire Year"
        if nodays1 <= 364:
            text = str(nodays1)+" days in each year"
        title = "Average Temperature (using TMID) for "+text+" between "+self.refstring
        plt.plot(self.all_years_mean[:,0],self.all_years_mean[:,1])
        plt.plot(trend_x,trend_y, color='black',  linestyle='dashed')
        
        plt.text(np.mean(trend_x),np.mean(trend_y)+1,str('Trend: '+str(round(slope*10,1))+"F/decade"),horizontalalignment='center')
        plt.title(title)
        plt.xlabel('Year')
        plt.ylabel('Temperature, Farenheit')
        #This creates a span for the baseline which is superimposed on the chart. 
        basespan_min=self.baseperiod[0,0]
        basespan_max=self.baseperiod[-1,0]
        plt.text((basespan_min + basespan_max)/2,np.nanmin(self.all_years_mean[:,1]),'Baseline',horizontalalignment='center')
        plt.axvspan(basespan_min, basespan_max, facecolor='g', alpha=0.5)
        #And the same for the reference period. 
        refspan_min=self.refperiod[0,0]
        refspan_max=self.refperiod[-1,0]
        if refspan_min == refspan_max :
            refspan_max = refspan_max + 1
        plt.axvspan(refspan_min, refspan_max, facecolor='r', alpha=0.5)
        plt.text((refspan_min + refspan_max)/2,np.nanmin(self.all_years_mean[:,1]),'Reference',horizontalalignment='center')
        
        plt.show()

        