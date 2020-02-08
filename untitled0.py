# -*- coding: utf-8 -*-
"""
Created on Sat Feb  8 13:12:27 2020

@author: abc73_000
"""
import io
a =[]
with open('C:\\Users\\abc73_000\\Desktop\\test.sam','rb') as f:
   
    while True:
         line = f.readline()
         if line:
             each= str(line.decode('ascii'))
             a.append(each)
         else:
            break
