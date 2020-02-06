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

#file_line = int(sys.argv[2]) 
with open ("C:/Users/abc73_000/Desktop/SRR7780979_DepthofCoverage_Distribution.txt",'r') as f:
    file = f.read()

table = file.split('\n')[:-1]

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

pos_number = []
freq_number =[]    
sum = 0    
line = 0
for con in table:
    content_list = makeOneSpace(con)
    freq = int(content_list[1])
    pos = int(content_list[2])
    pos_number.append(pos)
    freq_number.append(freq)
    total = freq*pos
    sum += total

df = {'Frequency': np.array(freq_number), 'Position': np.array(pos_number)}    
data_frame = pd.DataFrame(data=df )

plt.figure(figsize=(16,9))
sns.set(font_scale=3)
sns.distplot(freq_number, bins = len(freq_number),kde=False, axlabel= 'Frequency', color='orange')      
plt.ylim(0,100)
plt.xlim(0,1000)
plt.close()
#g.legend.remove() 
    
    
       
