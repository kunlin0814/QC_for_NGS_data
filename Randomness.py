# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 18:53:41 2020

The script take the depth of coverage as input file 
to identify the frequency of the coverage of depth for all of the base
position for each sample 


@author: abc73_000
"""

## look at all the positions to the end of the bases

import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import skew 
from scipy.stats import poisson
from sklearn.metrics import mean_squared_error
from math import sqrt



input_file = "/Users/kun-linho/Desktop/Pan_cancer_mapping_result/Distribution/Mammary/Normal/SRR7780741_DepthofCoverage_Distribution.txt"
#sys.argv[1]
#'G:\\Pan_cancer\\Pan_cancer_mapping_result\\Distribution\\Mammary\\Normal\\SRR7780741_DepthofCoverage_Distribution.txt'
#sys.argv[1]
file_name = "SRR7780741"
#sys.argv[2]



with open(input_file,'r')as f:
    file = f.read().split('\n')[:-1]

file.remove('1 Total_Depth')
summary = {}
original_list=[]
for i in range(len(file)):
    freq = int(file[i].split(" ")[0])
    pos = int(file[i].split(" ")[1])
    summary[pos]= int(freq)
    for j in range(freq):
        original_list.append(pos)

last_pos = int(list(summary.keys())[-1])

for i in range(last_pos+1):
    if i not in summary.keys():
        summary[i]=0
order_summary={}
for i in sorted(summary.keys()):
    order_summary[i]=summary[i]


total_data = pd.DataFrame(order_summary.items())
total_data.columns = [ 'Position', 'Frequency',]

freq_arr = total_data['Frequency'].values

average = np.mean(np.array(original_list)) 
mu = average

poisson_list=[ poisson.pmf(i,mu) for i in range(last_pos+1)]

    #for i in range(601):
#    poisson_list.append(poisson.pmf(i,mu))

rmse = sqrt(mean_squared_error(freq_arr, np.array(poisson_list)))    
sumOfSqerror= sqrt(sum((freq_arr-np.array(poisson_list))**2))

## need to check the longest continuous number


poisson_list_count=[poisson.pmf(i,mu)*Totalline for i in range(max_line)]

rmse_count = sqrt(mean_squared_error(freq_arr, np.array(poisson_list_count)))    
sumOferror_count= sqrt(sum((freq_arr-np.array(poisson_list_count))**2))
   
std = np.std(np.array(original_list))
   
#skewness = skew(np.array(original_list))

output = open(file_name+'_randomness_summary.txt','w')

output.write(file_name+'\t'+str(average)+'\t'+str(std)+'\t'+str(rmse)+'\t'+str(sumOferror)+'\t'+str(rmse_count)+'\t'+str(sumOferror_count)+'\n')
output.close()


## median= position (sum(freq)+1)/2 th position 
## Calcaulte Standard deviation ###
## sqrt(sum of (freq*(value - mean)**2)/sum of (freq)) 


"""
plt.figure(figsize=(16,9))
sns.set(font_scale=3)
sns.lineplot(x ='Position', y = 'Frequency', data =df)
#sns.distplot(freq_number, bins = len(freq_number),kde=False, axlabel= 'Frequency', color='orange')      

#plt.close() 

def makeOneSpace(String):
    content = String.split(' ')
    i = 0
    start = True
    while start:
        if content[1].isdigit() and content[2].isdigit():
            start = False
        else:
            i+=1
            content = String.split(' ')[i:]
            #print(i)   
    return content 
"""

