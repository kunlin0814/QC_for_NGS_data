# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 14:53:50 2020

WGS randomness will out of memory even use 40G,
So, I change a new version without created an original list
The result was compared with the GATK depth of coverage file "bed.sample_summary"
The mean coverage and the total read = total_pos_cover is the same
The whole thing will significantly speed out the process even for WGS


@author: abc73_000
"""
## look at all the positions to the end of the bases

import sys
import numpy as np
import pandas as pd
#import seaborn as sns
#import matplotlib.pyplot as plt
from scipy.stats import skew 
from scipy.stats import poisson
from sklearn.metrics import mean_squared_error
from math import sqrt
from collections import OrderedDict
#from matplotlib.backends.backend_pdf import PdfPages

"""
with open("/Users/kun-linho/Desktop/Pan_cancer_mapping_result/Distribution/Osteo/Normal/Normal_list")as f:
    file_list = f.read().split('\n')[:-1]
    
pp = PdfPages('/Users/kun-linho/Desktop/Normal_Osteosarcoma.pdf')
"""
input_file = sys.argv[1]
#"G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Validation_Data_Set\\MC\\004\\SRR9911376_DepthofCoverage_Distribution.txt"
file_name = sys.argv[2]
#'SRR10351814'
Cancer_type = sys.argv[3]
#'Mammary_Cancer'
Status = sys.argv[4]
#'Normal'

with open(input_file,'r')as f:
    file = f.read().split('\n')[:-1]
file.remove('1 Total_Depth')

last_pos = int(file[-1].split(" ")[1])
if last_pos >= 1000:
    last_pos = last_pos
else:
    last_pos = 1000 

toal_pos_cover=0
total_cover = 0
total_pos = 0
full_position_simultation ={}
for i in range(len(file)):
    number_position = int(file[i].split(" ")[0])
    cover_times = int(file[i].split(" ")[1])
    full_position_simultation[cover_times]= int(number_position)
    total_cover += cover_times
    total_pos +=number_position
    toal_pos_cover+=number_position*cover_times

ratio_order_dict=OrderedDict()
order_dict = OrderedDict()

## fill up the gap data, if there is no data at that position, then fill up with 0    
for i in range(last_pos+1):
    if i not in full_position_simultation.keys():
        ratio_order_dict[i]=0
        order_dict[i]=0
    else:
        ratio_order_dict[i] = full_position_simultation[i]/total_pos
        order_dict[i] = full_position_simultation[i]
            
mean_coverage = toal_pos_cover/total_pos
    

total_data = pd.DataFrame(list(order_dict.items()),columns=['number_position', 'cover_times'])

freq_arr = total_data['cover_times'].values
pos_arr = total_data['number_position'].values
mu = mean_coverage
prob_arr= freq_arr/total_pos

poisson_fract_list=[]

for i in range(last_pos+1):
    value = poisson.pmf(i,mu)
    poisson_fract_list.append(value)
    
rmse = sqrt(mean_squared_error(prob_arr, np.array(poisson_fract_list))) ## compare with the proportion of poisson dis

"""
poisson_list_count =[]
for i in range(last_pos+1):
    value = poisson.pmf(i,mu)*total_pos
    poisson_list_count.append(value)
"""

output = open(file_name+'_randomness_summary.txt','w')

output.write('file_name'+'\t'+'average'+ \
'\t'+'rmse'+'\t'+ \
'Cancer_type'+'\t'+ 'Status'+'\n')

output.write(file_name+'\t'+str(mean_coverage)+ \
'\t'+str(rmse)+'\t'+ \
+str(Cancer_type)+'\t' \
+str(Status)+'\n')

output.close()



"""
# alternative way to create an order an fill up 0
current_pos=0
i = 0
while i < len(file):
    number_position = int(file[i].split(" ")[0])
    cover_times = int(file[i].split(" ")[1])
    
    if current_pos == cover_times:
        full_position_simultation[cover_times] = number_position
        current_pos+=1
        i+=1    
    else:
        while current_pos < cover_times:
            full_position_simultation[current_pos] = 0
            current_pos+=1
        i+=1
        current_pos+=1
        full_position_simultation[cover_times]= number_position
        
        
# here is the plot section        
p = sns.lineplot(x ='Pos', y = 'Prob', data =Prob_data_1000)
plt.title(Cancer_type+"_"+file_name+"_"+Status)
p.set(xlabel='Coverage_Depth', ylabel= "fraction of the number of position")
#sns.distplot(total_data['Frequency'],hist=False, kde=True, axlabel= 'Frequency', color='orange')
p = sns.distplot(np.array(original_list),kde=False, axlabel= 'Position', color='black', bins=200)  
plt.title(file_name)
p.set_yscale('log')   # set into log scale 
#pp.savefig()
#p = p.map(plt.hist, "value", color="r", log=True)
plt.figure(figsize=(13,6))
sns.set(font_scale=2)      
#sns.kdeplot(np.array(original_list), shade=True,bw='scott');
#pp.savefig()
plt.close() 
#pp.close()
  
 
    
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

