# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 18:53:41 2020

The script take the depth of coverage as input file 
to identify the frequency of the coverage of depth for all of the base
position for each sample 

## look at all the positions to the end of the bases
## this script is used to create a plot that describes the sequence distriubtion.

## fill up the gap data, if there is no data at that position, then fill up with 0
@author: abc73_000
"""


import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import poisson
from sklearn.metrics import mean_squared_error
from math import sqrt
import glob, os
from matplotlib.backends.backend_pdf import PdfPages
from pathlib import Path
from math import sqrt
from collections import OrderedDict

pdf = PdfPages('C:/Users/abc73_000/Desktop/OSA.pdf')
location = "C:/Users/abc73_000/Desktop/Validation Data_Set/"
CancerType = "OSA"
#CaseName = "P1"
#file_folders = os.listdir(location+CancerType+"/"+"*"+"/"+"*txt")
for root, dirs, files in os.walk(location+CancerType+"/"):
    for file_name in files:
        input_file =  root+"/"+file_name
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
                
            
        total_data = pd.DataFrame(list(order_dict.items()),columns=['cover_times','number_position'])
        
        cover_times_arr = total_data['cover_times'].values
        number_pos_arr = total_data['number_position'].values
        prob_arr= number_pos_arr/total_pos
        
        cover_1000 = cover_times_arr[0:1000]
        pos_1000 = number_pos_arr[0:1000]
        ratio_1000 = prob_arr[0:1000]
        
        plot_input = {'cover_times':cover_1000,'number_position':pos_1000,'ratio_of_pos':ratio_1000}
        
        plot_data = pd.DataFrame(plot_input)
        
    
        plt.figure(figsize=(13,6))
        sns.set(font_scale=2)
        sns.lineplot(data = plot_data, x="cover_times",y="number_position")
        plt.title(file_name)
        #p.set_yscale('log')  # set into log scale 
        pdf.savefig()
        plt.close()
        
        plt.figure(figsize=(13,6))
        sns.set(font_scale=2)
        plt.title(file_name)      
        sns.lineplot(data = plot_data, x="cover_times",y="ratio_of_pos")
        
        pdf.savefig()
        plt.close()
        
        
pdf.close()
"""
for i in range(1000):
    value = poisson.pmf(i,mu)*pos_1000_line
    poisson_list_count_1000.append(value)
"""     


"""
output = open(file_name+'_randomness_summary.txt','w')

output.write(file_name+'\t'+str(average)+'\t' \
+str(std)+'\t'+str(rmse)+'\t'+str(sumOfSqerror)+'\t' \
+str(Cancer_type)+'\t'+str(Status)+'\n')

output.close()
"""

"""
p = sns.lineplot(x ='Pos', y = 'Prob', data =Prob_data_1000)
plt.title(Cancer_type+"_"+file_name+"_"+Status)
p.set(xlabel='Coverage_Depth', ylabel= "fraction of the number of position")
sns.distplot(total_data['Frequency'],hist=False, kde=False, axlabel= 'Frequency', color='orange')
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
   
  
    
plt.figure(figsize=(13,6))
sns.set(font_scale=2)
#p = sns.lineplot(x ='Position', y = 'Frequency', data =total_data)
sns.distplot(total_data['Frequency'],hist=False, kde=True, axlabel= 'Frequency', color='orange')
p = sns.distplot(np.array(original_list),kde=False, axlabel= 'Position', color='black', bins=200)  
plt.title(file_name)
p.set_yscale('log')   # set into log scale 
pp.savefig()
#p = p.map(plt.hist, "value", color="r", log=True)
plt.figure(figsize=(13,6))
sns.set(font_scale=2)      
sns.kdeplot(np.array(original_list), shade=True,bw='scott');
pp.savefig()
plt.close() 
pp.close()
"""



