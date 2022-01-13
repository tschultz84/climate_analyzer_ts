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

dataval=bzcalc.all_years_mean

#This function creates a set of slopes by permuting the data.
#Data must be in the format of self.all_years_mean,
#so that it can be submitted to self.create_trend_line
#n is the number of permutations. 
#If hist is true, show a histogram.
#It returns the precentile value of the real slope in the permuted set.

def permuted_slopes(data,n=1000,hist=False,display=False):
    if display == True : startTime = time.time()
    #First, create the actual values. 
    real_data = bzcalc.create_trend_line(data,array=False)
    #initialize the array that wil be returned. 
    slopes=np.array([real_data[0]])
    #Create a loop over n, to create the permutations.
    for i in np.arange(0,n):
        #Pull out the data into seaprate arrays. 
        years = data[:,0]
        values=data[:,1]
        #Permute the values data. 
        values = np.random.permutation(values)
        #Combine them into a format suitable for submissin to create_trend_line
        comb = np.array([years,values]).reshape(-1,2)
        next1 = bzcalc.create_trend_line(comb,array=False)
        #Then append. You multiply next1 * 10 in orer toc reate
        #an output in the units of F / decade.
        slopes = np.append(slopes,10*next1[0])
    if hist == True:
        #Creates a histogram of all permuted slope values. 
        plt.hist(slopes,bins=30,density=True)
        plt.title('Frequency of permuted Slopes (line shows real slope')
        plt.ylabel('Fraction of Total')
        plt.xlabel('Slope Value (F warming per decade')
        plt.axvline (x=real_data[0],color='k', linestyle='--')
        plt.show()
    
    #This calculates the percentile value. ()
    #This is the actual value that indicates the probability of this happening 
    #by chance. 
    percentile_value = np.count_nonzero(slopes<real_data[0])/len(slopes)   

    if(display):
        executionTime = (time.time() - startTime)
        print("Time taken to calculate "+str(n)+" permuted slope values." + str(executionTime))   
        print("That is "+str(round(100*executionTime/n,2))+" seconds per 100 permutations.")
    return percentile_value
orig = dataval[:,1]
test=permuted_slopes(dataval,5000,True,True)
print(test)
