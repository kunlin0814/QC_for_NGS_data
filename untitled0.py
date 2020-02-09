# -*- coding: utf-8 -*-
"""
Created on Sat Feb  8 13:12:27 2020

@author: abc73_000
"""
import io
from datetime import datetime
import numpy as np

now = datetime. now()

sam_file="C:\\Users\\abc73_000\\Desktop\\CMT-SRR7780976-test2.sam"
#sys.argv[1]
file_name="CMT-SRR7780976"
#sys.argv[2]
read_length = 101
#int(sys.argv[3])



## criteria is Reads start site <=Exon end site and exon start site <= reads end site
def withinRegion(reads_position, exom_start, exom_end):
     reads_end = reads_position + read_length-1
     reads_start = reads_position
     if (reads_start <= exom_end and exom_start <= reads_end):
         return 1
     elif (reads_start > exom_end):
         return 0
     elif (reads_end < exom_start):
         return -1


def binarySearch (arr, left, right, reads_position, read_length): 
    
    # Check base case 
    if right >= left: 
        mid = int((left+right)/2)
        exon_start_value = int (arr[mid][0])
        exon_end_value = int (arr[mid][1])
        # If element is present at the middle itself 
        if (withinRegion(reads_position, exon_start_value, exon_end_value ) == 1):
            return mid 
          
        # If element is smaller than mid, then it can only 
        # be present in left subarray 
        elif (withinRegion(reads_position, exon_start_value, exon_end_value ) == -1):
            return binarySearch(arr, left, mid-1, reads_position, read_length) 
  
        # Else the element can only be present in right subarray 
        elif (withinRegion(reads_position, exon_start_value, exon_end_value ) == 0) : 
            return binarySearch(arr, mid+1, right, reads_position, read_length) 
  
    else: 
        # Element is not present in the array 
        return -2


#small_interval_dict ={}
#large_interval_dict ={}
Total_interval_dict = {}
with open('G:\\Pan_cancer\\Mapping_source\\Canis_familiaris.CanFam3.1.81.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list', 'r') as f:
    file = f.read()


CDS = file.split('\n')[:-1]
    
for i in range(len(CDS)):
    chrom = CDS[i].split(':')[0]
    start = int(CDS[i].split(':')[1].split('-')[0])
    end = int(CDS[i].split(':')[1].split('-')[1])
    if chrom not in Total_interval_dict.keys():
        Total_interval_dict[chrom]=[(start, end)]
    else:
        Total_interval_dict[chrom].append((start, end))
   

for i in Total_interval_dict.keys():
    Total_interval_dict[i]=np.array(Total_interval_dict[i],dtype=np.int32)
    np.sort(Total_interval_dict[i])
   

unique = 0 #3
#duplicate = 0 #3
#Onemapped = 0 #5,9
#incorrect = 0 #1
#unmapped = 0 #13

transcript_list =[]    
total = 0 
pass_line =0


with open(sam_file,'rb') as f:   
     while True:
         line = f.readline()
         if line:
            
            file_lst=line.decode("utf-8").split('\t')
            
            if '@' in file_lst[0]:
                pass
            else :
                total += 1
                reads_name = file_lst[0]
                reads_chr = file_lst[2]
                reads_position = int(file_lst[3])
                status = int(file_lst[1])%16
                if status == 3:
                    for ele in file_lst:
                        if 'XT:' in ele:
                            status2 = ele.split(':')[2]
                            if status2 == 'U' or status2 == 'M':
                                unique += 1
                                exome_loc = Total_interval_dict[reads_chr]
                                if (binarySearch(exome_loc,0, len(exome_loc)-1,reads_position,read_length)!=-2):
                                    transcript_list.append(reads_name)
                                    pass_line+=1
         else:
            break
                                
                           

#dup = [key for (key, value) in Counter(transcript_list).items() if value > 1 and key]
pairs = total / 2
current_time = datetime.now()
print(current_time-now)


e =[None]*1000000
now =datetime.now()
for i in range(100000):
    e[i]= 'SRR100'+str(i)
    
a = np.array(e,dtype=np.unicode_)   
current_time = datetime.now()
print(current_time-now)
        

now =datetime.now()
b= ['SRR100'+str(i) for i in range(100000)]
current_time = datetime.now()
print(current_time-now)
         




"""
a =[]
with open(sam_file,'rb') as f:
  for i in f:
      a.append(i)
"""