#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 17:01:36 2020

@author: kun-linho
"""

import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
#from scipy.stats import skew 
#from scipy.stats import poisson
#from sklearn.metrics import mean_squared_error
from math import sqrt

input_file = sys.argv[1]
file_name = sys.argv[2]


mapping = pd.read_csv(file_name+"_mapping_quality")
#ile_name+"_mapping_quality")
value = mapping.iloc[0:].values
mean = np.mean(value)
std = np.std(value)
median = np.median(value)

result = open(file_name+'_mapping_quality_summary.txt','w')

result.write(str(mean)+'\t'+str(median)+'\t'+str(std)+'\n')
result.close()

"""
plt.figure(figsize=(16,9),dpi=300)
sns.set(font_scale=1)
plt.xlim(0,80)
sns.distplot(mapping,kde=False)
plt.savefig('/Users/kun-linho/Desktop/SRR7781089_mapping_quality.png',format='png',dpi=300)
#sns.distplot(freq_number, bins = len(freq_number),kde=False, axlabel= 'Frequency', color='orange') 
plt.close()
"""