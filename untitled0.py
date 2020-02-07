#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 09:40:32 2020

@author: kun-linho
"""

def main():
    car ={'plate':"GA", 'millage':50000}
    car['brand']='WA'
    car['year']=2015
    car['owner']='lauren'
    car['seller']= 'VW Athens'
    print(car)
    for i,j in car.items():
        print(i,j )
    for i in car.keys():
        print (i)
    for i in car:
        print (car[i])
    print(len(car))
    car.pop('seller')
    print(len(car),car)
main()
        