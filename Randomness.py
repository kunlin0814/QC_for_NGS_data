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

Total_line = int(sys.argv[2]) 
with open ("/Users/kun-linho/Desktop/test_distribution",'r') as f:
    file = f.read()

table = file.split('\n')[:-1]
table.remove('1 Total_Depth')

   

pos_number = []
freq_number =[]    
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

average = sum/Total_line

Std = 


### Calcaulte Standard deviation ###
## sqrt(sum of (freq*(value - mean)**2)/sum of (freq)) 

df = {'Frequency': np.log2(np.array(freq_number)), 'Position': np.array(pos_number)}    
data_frame = pd.DataFrame(data=df )

plt.figure(figsize=(16,9))
sns.set(font_scale=3)
sns.lineplot(x ='Position', y = 'Frequency', data =df)
#sns.distplot(freq_number, bins = len(freq_number),kde=False, axlabel= 'Frequency', color='orange')      

#plt.close()
#g.legend.remove() 
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

