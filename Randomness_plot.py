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
from scipy.stats import skew 
from scipy.stats import poisson
from sklearn.metrics import mean_squared_error
from math import sqrt
from matplotlib.backends.backend_pdf import PdfPages
import glob, os
"""
with open("/Users/kun-linho/Desktop/Pan_cancer_mapping_result/Distribution/Osteo/Normal/Normal_list")as f:
    file_list = f.read().split('\n')[:-1]
    

"""

pp = PdfPages('C:/Users/abc73_000/Desktop/OSA.pdf')
location = "C:/Users/abc73_000/Desktop/Validation Data_Set/"
CancerType = "OSA"
#CaseName = "P1"
#file_folders = os.listdir(location+CancerType+"/"+"*"+"/"+"*txt")
for root, dirs, files in os.walk(location+CancerType+"/"):
    for file_name in files:
        input_file =  root+"/"+file_name
        #print(input_file)
        
        #file_name = "SRR9911376_DepthofCoverage_Distribution.txt"
        #sys.argv[1]
        
        #i.split('_')[0]
        #sys.argv[2]
        #"SRR7780741"
        #sys.argv[2]
        #'Mammary_Cancer'
        #sys.argv[3]
        Status =  'Normal'
        #'Normal'
        #sys.argv[4]
        ## fill up the gap data, if there is no data at that position, then fill up with 0
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
                
        total_line = len(original_list)
        
        if original_list[-1] > 1000:
            last_pos = original_list[-1]
        else:
            last_pos = 1000
            
        
        
        for i in range(last_pos+1):
            if i not in summary.keys():
                summary[i]=0
                
        order_summary={}
        for i in sorted(summary.keys()):
            order_summary[i]=summary[i]
        
        total_data = pd.DataFrame(list(order_summary.items()),columns=['Position', 'Frequency'])
        
        freq_arr = total_data['Frequency'].values
        pos_arr = total_data['Position'].values
        freq_arr_1000 = freq_arr[0:1000]
        pos_arr_1000 = pos_arr[0:1000]
       
        total_array = np.array(original_list)
        
        
        original_list_1000=[]
        for i in range(1000):
            pos = pos_arr[i]
            freq = freq_arr_1000[i]
            for j in range(freq):
                original_list_1000.append(pos)
        
        
        #average = np.mean(np.array(original_list))
       
        #std = np.std(np.array(original_list))
       
        
        #mu = average
        pos_1000_line = len(original_list_1000)
        #prob_arr= freq_arr/total_line
        prob_arr_1000 = freq_arr_1000/ pos_1000_line
        #freq_list= list(freq_arr))
        
        Prob_data_1000 = pd.DataFrame(data= {'Prob' :prob_arr_1000,'Pos': pos_arr_1000})
        pd.DataFrame(list(order_summary.items()),columns=['Position', 'Frequency'])
        #rmse = sqrt(mean_squared_error(prob_arr, np.array(poisson_fract_list))) ## compare the proportion
        #sumOfSqerror= sqrt(sum((prob_arr-np.array(poisson_fract_list))**2)) ## rmse that didn't divide the N, vs proportion
 
        
        poisson_list_count =[]
        poisson_list_count_1000 =[]
        
        
        
        
        plt.figure(figsize=(13,6))
        sns.set(font_scale=2)
        p = sns.distplot(np.array(original_list_1000), hist=True,kde=False, \
                     color = 'darkblue', \
                     bins = 500, \
                     hist_kws={'edgecolor':'black'}, \
                     kde_kws={'linewidth':3});
                     #kde_kws={'linewidth': 3})
        plt.title(file_name)
        p.set_yscale('log')   # set into log scale 
        
        pp.savefig()
        plt.close()
        
        plt.figure(figsize=(13,6))
        sns.set(font_scale=2)
        plt.title(file_name)      
        sns.kdeplot(np.array(original_list_1000), shade=True,bw=.1);
        
        pp.savefig()
        plt.close()
        
pp.close()
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



