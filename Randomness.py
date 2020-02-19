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


input_file = "/Users/kun-linho/Desktop/Pan_cancer_mapping_result/Distribution/Mammary/Normal/SRR7780741_DepthofCoverage_Distribution.txt"
#sys.argv[1]
#'G:\\Pan_cancer\\Pan_cancer_mapping_result\\Distribution\\Mammary\\Normal\\SRR7780741_DepthofCoverage_Distribution.txt'
#sys.argv[1]
file_name = "SRR7780741"
#sys.argv[2]

#sys.argv[1]
#'G:\\Pan_cancer\\Pan_cancer_mapping_result\\Distribution\\Mammary\\Normal\\SRR7780741_DepthofCoverage_Distribution.txt'
#sys.argv[1]
#sys.argv[2]
#'SRR7780741'
#sys.argv[2]
with open (input_file,'r') as f:
    file = f.read()

table = file.split('\n')[:-1]
table.remove('1 Total_Depth')
   

pos_number = []
freq_number =[]
original_list=[]    
Total_line = 0    
#line = 0
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

total_data = pd.read_csv("/Users/kun-linho/Desktop/Pan_cancer_mapping_result/Distribution/Mammary/Normal/SRR7780741_DepthofCoverage_Distribution.txt", sep=" ", header = None)
final_total = total_data.drop(1,axis = 0)



#frequency_array = np.array(freq_number)/Total_line
df = {'Frequency': np.array(freq_number), 'Position': np.array(pos_number)}    
data_frame = pd.DataFrame(data=df )


last_Position = data_frame['Position'].iloc[-1]

#poisson_number = int(data_frame.iloc[-1]["Position"])
Total600line = sum(data_frame['Frequency'][0:600]*data_frame['Position'][0:600])
event_arr = (data_frame['Frequency'][0:600].values)/Total600line
freq_arr = data_frame['Frequency'].values[0:600]

poisson_number = data_frame["Position"].values

max_line= 600
        

Totalline = sum(data_frame['Frequency'][0:max_line]*data_frame['Position'][0:max_line])
event_arr = (data_frame['Frequency'][0:max_line].values)/Totalline
freq_arr = data_frame['Frequency'].values[0:max_line]

average = np.mean(np.array(original_list)) 
mu = average

poisson_list=[ poisson.pmf(i,mu) for i in range(max_line)]

    #for i in range(601):
#    poisson_list.append(poisson.pmf(i,mu))

rmse = sqrt(mean_squared_error(event_arr, np.array(poisson_list)))    
sumOferror= sqrt(sum((event_arr-np.array(poisson_list))**2))

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

