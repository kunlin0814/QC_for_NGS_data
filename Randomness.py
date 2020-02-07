# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 18:53:41 2020

@author: abc73_000
"""

import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import skew 


input_file = sys.argv[1]
file_name = sys.argv[2]
with open (input_file,'r') as f:
    file = f.read()

table = file.split('\n')[:-1]
table.remove('1 Total_Depth')
   

pos_number = []
freq_number =[]
original_list=[]    
sum = 0    
line = 0
for content in table:
    #content_list = makeOneSpace(con)
    content_list = content.split()
    freq = int(content_list[0])
    pos = int(content_list[1])
    pos_number.append(pos)
    freq_number.append(freq)
    total = freq*pos
    sum += total
    for i in range(freq):
        original_list.append(pos)
          
std = np.std(np.array(original_list))
average = np.mean(np.array(original_list))    
skewness = skew(np.array(original_list))

output = open(file_name+'_randomness_summary.txt','w')

output.write(file_name+'\t'+str(average)+'\t'+str(std)+'\t'+str(skewness)+'\n')
output.close()


## median= position (sum(freq)+1)/2 th position 
## Calcaulte Standard deviation ###
## sqrt(sum of (freq*(value - mean)**2)/sum of (freq)) 
"""
df = {'Frequency': np.log2(np.array(freq_number)), 'Position': np.array(pos_number)}    
data_frame = pd.DataFrame(data=df )

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

