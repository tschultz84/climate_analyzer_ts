# -*- coding: utf-8 -*-

import pandas as pd
import sys
import yaml
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
from sklearn.linear_model import LinearRegression
from scipy.stats.stats import pearsonr
from IPython.display import display_html 


#All analysis fucntions for a station are done in this object.
#The input is stationdata, which shoudl be self.station_data from the StationLoad object
#the refst and refend are both strings, marking the beginning and end of the reference period of analysis

#display prints the function outputs

class StationAnalyzer :
    def __init__(self,stationdata,refst='2020-01-31',refend='2020-12-31',display=False):
        
        self.display=display
        #This YAML file contains a great deal of static information, 
        #such as directory information. 
       # yaml_dir="C:\\Users\\14154\\OneDrive\\Python\\climate_mapper\\python\\climate_analyzer_ts\\"
        yaml_file = open("load_stats_static.yaml")
        self.yaml = yaml.load(yaml_file, Loader=yaml.FullLoader)
    
        #Creating the assigned value containing every single element. 
        self.all_data_long=stationdata
        #Then assigning off TMID, TMIN, and TMAX values, which are different
        #columns of the input data.
        self.tmid_array = self.all_data_long[:,[0,1,2,3,6]]
        self.tmin_array = self.all_data_long[:,[0,1,2,3,5]]
        self.tmax_array = self.all_data_long[:,[0,1,2,3,4]]
        
        #This then creates the yearly averages.
        #First, find the beginning and end years, and the number of years total.
        self.alltime_startyear=np.nanmin(self.tmax_array[:,0])
        self.alltime_endyear=np.nanmax(self.tmax_array[:,0])
        self.alltime_number_years=self.alltime_endyear-self.alltime_startyear
        
        #Creates strings of the ref periods for later access. 
        self.refststr=refst
        self.refendst=refend
        
        """This creates an array of annual average temperatures across the entire
        time period.
        """
        yearsout=np.empty([0,2])
        uniqueyears=np.unique(self.tmid_array[:,0])
        for year in uniqueyears[:-1]:   
            #Select only the index values wehre the year is equal to year.
            index=np.where(self.tmid_array[:,0]==year)
            #Then, selects data for this year
            #All columns are pulled, since we have to filter to ensure ther eare 
            #12 months for every year evaluated. 
            all_years_data=self.tmid_array[index]
            
            #And then evaluates their mean. 
            all_years_means=np.nanmean(all_years_data[:,4])
            row=np.array([year,all_years_means])
            yearsout=np.append(yearsout,[row],axis=0)
        #Then assigns a variable to the output array.     
        self.all_years_mean=yearsout
        
        """These calls create arrays, containing a number of entries in the form
        [year,day of year] which indicates all days in the references and baseline
        periods. 
        The TMID array values are then defined in these periods.        """
        
        #Creating the arrayed reference period. 
        self.refperiod=self.find_ref_range(refst,refend)
        #creates baseline range.
        self.find_baseline_range()
        
        #Creatse the actual values of TMID in the periods.
        self.tmid_ref_data=self.tmix_selection(self.refperiod,self.tmid_array)
        self.tmid_base_data=self.tmix_selection(self.baseperiod,self.tmid_array)
        
        """#This creates an array with monthyl average temperatures. """
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
        self.all_months_mean=monthout
        
        #Runs the key_metrics function if desired.
        if self.display==True:
            #Generate variables. 
            self.key_metrics()
            #show everything in a nice format. 
            display_html(self.key_metrics_table)
            display_html(self.key_stats)
            self.key_charts()
            
            
    """BEGIN FUNCTIONS"""    
    
    """DATA CLEANING AND WRANGLING FUNCTIONS
    The following set of functions make various transformations on data
    and create new arrays in order to complete various analyses.
    """  
        
    #This function finds the Day of year, assuming every month has 31 days for simplicity
    #Its input is a pd.DateTime object, it outputs an integer.
    def find_doy(self,pddate):
        doy = pddate.day+(pddate.month-1)*31
        return doy[0]
    
    #This function creates beginning and end dates for two different pd datetimeIndex
    #objects, and creates a range of days, int he form of [[Year1,Day of Year1],[Year1,DOY2],...]
    #It takes in two pd.DatetimeIndex objects which are the beginning and end of the reference
    #periods.
    #It is used to create the reference period range of days.
    def find_ref_range(self,pddate1,pddate2):
        #This function is supposed to be passed a Pandas datetimeIndex.
        #If it's a string, a conversion is made. 
        if type(pddate1)==str:
            pddate1=pd.DatetimeIndex([pddate1])
        if type(pddate2)==str:
            pddate2=pd.DatetimeIndex([pddate2])
        
        #First, create simple forms of the Day of year and Year. 
        doy1 = self.find_doy(pddate1)
        doy2 = self.find_doy(pddate2)
        year1=pddate1.year[0]
        year2=pddate2.year[0]
        #This simply throws an error if the years are not ordered properly.   
        if year1>year2:        
          print("Your reference period is improperly defined.")
          print("THe first year must be before the second year.")
          sys.exit("Break Error due to bad dates. .")
        #First, define if the reference period beginning and end are in the same eyar. 
        if year1==year2:
            #Raise an exception if the reference period end is before the beginning.
            if doy1>doy2:
                
                print("Your reference period is improperly defined.")
                print("THe first date must be before the second date.")
                sys.exit("Break Error due to bad dates. .")
            #But if the beginning and end are in the same year, then the array is very simple, 
            #jsut each row corresonding to the same year,
            #with a year and Day of Year.
            alldays = np.arange(doy1,doy2+1)
            #allyears = year1
            returner = np.transpose(np.array([np.repeat(year1,len(alldays)),alldays]))
            
        #If year1 is before year2, the array is filled with day of year valuse
        #that include all days in the year. 
        if year1 < year2:
            #This prints a message, indicating that reference periods including mroe than one year
            #include the entire period of time.
            print("You included more than one year in the reference period.")
            print("The analysis will include all days in each year included.")
            print("In order to avoid seasonal bias, I have to average over the entire year.")
            print("For example if you entered 2020 Jan 1 to 2021 July 1 as a reference,")
            print("this biases your result to be much warmer than average years")
           
           #Define the returned array, which is empty.
            returner = np.empty([0,2])
            #First, define all the days in the list, which includes all days in the year
            dayskarr = np.arange(1,372+1)
            #for k in np.arange(year1+1,year2):
            for k in np.arange(year1,year2+1):
                #Then, over all years in between year 1 and year 2, create a new array
                #This one has DOY from 1 to 372 and every year in the intermediate period.               
                yearkarr = np.transpose(np.array([np.repeat(k,len(dayskarr)),dayskarr]))
                #This then appends to returner this array.
                returner = np.append(returner,yearkarr,axis=0)    
       
        return returner    
    
    #This then finds the Days of Year for the analyzed baseline period, using
    #the same days of year as the reference period, to ensure the comparison is fair.
    def find_baseline_range(self):
        #Takes the unique days from the refperiod.
        baseline_doy = np.unique(self.refperiod[:,1])
        #Then, all the baseline years.
        baseline_years = np.unique(self.all_years_mean[:,0])
        #and finally, limits it to the first set of baseline years.
        baseline_years=baseline_years[0:self.yaml['BASENOYEARS']]
        
        #Then, create an array which contains all of the elemnts. 
        #THe first column are the years.
        #The second the DOY. 
        #Creates a template. 
        returner = np.empty([0,2])
        #Loops over all unique years. 
        for k in np.arange(np.min(baseline_years),np.max(baseline_years)):
            #THen, creates each row element and appends it.             
            yearkarr = np.transpose(np.array([np.repeat(k,len(baseline_doy)),baseline_doy]))
            returner = np.append(returner,yearkarr,axis=0)
        #This creates a formatted string defined the range of the baseline period.
        
        #The following lines turn the baseline period provided into a nicely formatted
        #string, for display in the final returner table.
        #First, find the first years in the base period, and make them into a string.
        firstbaseyear = returner[0,0]
        lastbaseyear = returner[-1,0]
        basepdyears = str(int(firstbaseyear))+" to "+str(int(lastbaseyear))
        
        #Then, select the DOY in year 1.
        basefirstyeardoy = returner[np.where(returner[:,0]==firstbaseyear)][:,1]
        #Extract the first and last DOY.
        firstbasedoy = basefirstyeardoy[0]
        lastbasedoy = basefirstyeardoy[-1]
        
        #find exemplars at this DOY.
        rowex1 = self.tmid_array[np.where(self.tmid_array[:,3]==firstbasedoy)][0]
        rowex2 = self.tmid_array[np.where(self.tmid_array[:,3]==lastbasedoy)][0]
        
        #Find the month and day for each.
        firstbasemo = rowex1[1]
        firstbaseday= rowex1[2]
        
        lastbasemo = rowex2[1]
        lastbaseday = rowex2[2]
        
        #Then create a string.
        modaystr = str(int(firstbasemo))+"-"+str(int(firstbaseday))+" to "+str(int(lastbasemo))+"-"+str(int(lastbaseday))
        #This creates an object out of the string for later reference.
        self.base_period_string = modaystr+" in "+basepdyears    
        #this is the baseline period used in calculations.
        self.baseperiod = returner
        
    
    #This function selects the days corresponding to "range1" from the dataset "data".
    #range1 is an array in the format of self.refperiod, a long list of rows in the form [year,doy]
    #data is  data, in the format of self.tmid_array
    def tmix_selection(self,range1,data):
        #Rename the input variables, which are subject to change.
        out=data
        ref=range1
        
        #Select only the years.
        year=out[:,0]

        #Select only years in the data which are present in range1.
        index = np.in1d(year,ref[:,0])
        #Select only the subset of years defined by index.
        out=out[index]
        #Select only the days of year in ref.
        index=np.in1d(out[:,3],ref[:,1])
        out=out[index]
        return out
    
     #This takes a row from tmid_array, tmax_array, etc., and returns a formatted string.
     #The excludeyear field excludes the year, which makes displaying the baseline period easier.
    def date_from_row(self,row,excludeyear=False): 
        year = str(int(row[0]))  
        month = str(int(row[1]))
        day = str(int(row[2]))
        if excludeyear==False:
            return year+"-"+month+"-"+day  
        if excludeyear==True:
            return month+"-"+day  
        
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
        recentyears=self.yaml['RECENT_TREND_YEARS']
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
        alltimestr1 = self.date_from_row(self.tmid_array[0])
        alltimestr2 = self.date_from_row(self.tmid_array[-1])
        #Creates strings describing the range of dates.
        refstr1 = self.date_from_row(self.tmid_ref_data[0])
        refstr2 = self.date_from_row(self.tmid_ref_data[-1])
        #This creates a dataframe of data points.
        returner = pd.DataFrame(
            {
             "Period of Measure":
                 ["All Time:" +alltimestr1+" to "+alltimestr2,
                 "Reference Period: "+refstr1+" to "+refstr2,
                 "Baseline Period: "+self.base_period_string],
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
        #print(" with a probability that this occured by chance less than "+str(self.yaml['ALPHA']))
        #print(" i.e., alpha is "+str(self.yaml['ALPHA'])+str("?"))
     
        if pvalue >= self.yaml['ALPHA']:
           #print("No.")
           self.ref_stat_sig = False
        if pvalue < self.yaml['ALPHA']:
           #print("Yes.") 
           self.ref_stat_sig = True

        trend_data = self.create_trend_line(self.all_years_mean,False)
        #print("The warming trend over the past "+str(self.yaml['RECENT_TREND_YEARS'])+" is: "+str(trend_data[0]*10)+" F degrees per decade.")
        #This creates a variable to pull out the warming, in F per decade.
        self.trend_data_str=str(round(trend_data[0][0]*10,2))+" Farenheit per decade"
        #This calculates the correlation coefficient and accompanying p-value
        #which states whether this correlation is statistically significant. 
        
        #First, format, and take only the recent years of data.
        recentyears=self.yaml['RECENT_TREND_YEARS']
        x5 = self.all_years_mean[-recentyears:,0]
        y5 = self.all_years_mean[-recentyears:,1]
        
        pearsond1 = pearsonr(x5,y5)
                
       # print("Outcome of Pearson Correlation Coefficient Analysis: ")
       # print("Coefficient of Correlation (R): "+str(pearsond1[0]))
       # print("Accompanying P value: "+str(100*pearsond1[1]))
        if pearsond1[1] < self.yaml['ALPHA']:
            is_real = True
        if pearsond1[1] >= self.yaml['ALPHA']:
            is_real = False
       # print("Is this trened real using Pearson correlation analysis? "+str(is_real))
        #Variables to pull out.
        self.trend_p = pearsond1[1]
        self.trend_is_real = is_real
        #This creates a dataframe of statistical information points.
        key_stats = pd.DataFrame(
            {
                "Metric":
                    ["Reference Minus Basline Temperature Change",
                     str(self.yaml['RECENT_TREND_YEARS'])+ " year warming trend"
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
                         [self.yaml['ALPHA'],
                          self.yaml['ALPHA']]
                
                }
        )
        #print(key_stats)
        
        self.key_metrics_table =returner
        self.key_stats = key_stats
        return returner   
    
    #end key_metrics
              
     #This creates several charts of interest.
    def key_charts(self):
  
        
        #Creates a histogram of daily TMID values. 
        ref_hist=self.tmid_ref_data[:,4]
        base_hist=self.tmid_base_data[:,4]
        plt.hist(ref_hist,bins=30,density=True,label="Reference Period", alpha=0.5)
        plt.hist(base_hist,bins=30,density=True,label="Base Period", alpha=0.5)
        plt.title('Freqency Histogram of TMID over Reference and Base Period')
        plt.xlabel('Fraction of Total at this TMID')
        plt.ylabel('Temperature, F (TMID)')
        plt.legend()
        plt.show()
        
            
        #Creates a plot over one year of the monthyl average TMID temperatures. 
        plt.plot(self.all_months_mean[:,0],self.all_months_mean[:,1])
        plt.title('Monthly Average of TMID')
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
        
        #Creates a plot over time of the annual average values. 
        plt.plot(self.all_years_mean[:,0],self.all_years_mean[:,1])
        plt.plot(trend_x,trend_y, color='black',  linestyle='dashed')
        plt.text(np.mean(trend_x),np.mean(trend_y)+1,str('Trend: '+str(round(slope*10,1))+"F/decade"),horizontalalignment='center')
        plt.title('Annual Average of Monthly TMID')
        plt.xlabel('Year')
        plt.ylabel('Annual Average Temperature, F (TMID)')
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

        