# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 18:53:41 2020

The script take the depth of coverage as input file 
to identify the frequency of the coverage of depth for all of the base
position for each sample 


@author: abc73_000
"""

## only look at the position to 600 base, account for 98.6% lines of total_lines
import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import skew 
from scipy.stats import poisson
from sklearn.metrics import mean_squared_error
from math import sqrt

input_file = '/Users/kun-linho/Desktop/Pan_cancer_mapping_result/Distribution/Mammary/Normal/SRR7780741_DepthofCoverage_Distribution.txt'
#sys.argv[1]
file_name = 'SRR7780741'
#sys.argv[2]
with open (input_file,'r') as f:
    file = f.read()

table = file.split('\n')[:-1]
table.remove('1 Total_Depth')
   

pos_number = []
freq_number =[]
original_list=[]    
Total_line = 0    
line = 0
for content in table:
    #content_list = makeOneSpace(con)
    content_list = content.split()
    freq = int(content_list[0])
    pos = int(content_list[1])
    pos_number.append(pos)
    freq_number.append(freq)
    total = freq*pos
    Total_line += total
    for i in range(freq):
        original_list.append(pos)

total_pos_arr = np.array(original_list)
original_list.clear()

#frequency_array = np.array(freq_number)/Total_line
df = {'Frequency': np.array(freq_number), 'Position': np.array(pos_number)}    
data_frame = pd.DataFrame(data=df )
#poisson_number = int(data_frame.iloc[-1]["Position"])
Total600line = sum(data_frame['Frequency'][0:601]*data_frame['Position'][0:601])
event_arr = data_frame['Frequency'][0:601].values/Total600line
mu = 100
poisson_list=[]
for i in range(601):
    poisson_list.append(poisson.pmf(i,mu))

rmse = sqrt(mean_squared_error(event_arr, np.array(poisson_list)))    
sumOferror= sum(event_arr[0:601]-np.array(poisson_list)[0:601])**2
   
std = np.std(np.array(original_list))
average = np.mean(np.array(original_list))    
skewness = skew(np.array(original_list))

output = open(file_name+'_randomness_summary.txt','w')

output.write(file_name+'\t'+str(average)+'\t'+str(std)+'\t'+str(skewness)+'\n')
output.close()


## median= position (sum(freq)+1)/2 th position 
## Calcaulte Standard deviation ###
## sqrt(sum of (freq*(value - mean)**2)/sum of (freq)) 



plt.figure(figsize=(16,9))
sns.set(font_scale=3)
sns.lineplot(x ='Position', y = 'Frequency', data =df)
#sns.distplot(freq_number, bins = len(freq_number),kde=False, axlabel= 'Frequency', color='orange')      

#plt.close() 
""" 
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

