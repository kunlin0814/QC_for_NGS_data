# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 11:16:45 2020

@author: abc73_000
"""
<<<<<<< HEAD
#### Windows system, please don't use multiprocessing ###

=======
#### windows system, don't use multiprocessing ###


### simplified verion
import multiprocessing as mp

def process_wrapper(lineByte):
    final_status =[]
    with open("/Volumes/Research_Data/MAC_Research_Data/Pan_cancer/test.sam") as f:
        f.seek(lineByte) ## seek(0) will move to the begining line
        line = f.readline().split('\t')
       
        
        if '@' in line[0]:
            pass
        else :
            status = int(line[1])%16
            final_status.append(status)
    return final_status 
        
#init objects
pool = mp.Pool(4)
jobs = []

#create jobs
if __name__ == "__main__":
    with open("/Volumes/Research_Data/MAC_Research_Data/Pan_cancer/test.sam") as f:
         ## tell you what is the position of the line to start to be processed
        for line in f:
        #pool.map(process_wrapper, f)
            jobs.append( pool.apply_async(process_wrapper, 0) )
            #nextLineByte = f.tell()
            #print(nextLineByte)
        
    #wait for all jobs to finish
    for job in jobs:
        job.get()
    
    #clean up
    pool.close()




##
import multiprocessing as mp,os

def process_wrapper(chunkStart, chunkSize):
    with open("input.txt") as f:
        f.seek(chunkStart)
        lines = f.read(chunkSize).splitlines()
        for line in lines:
            process(line)

def chunkify(fname,size=1024*1024):
    fileEnd = os.path.getsize(fname)
    with open(fname,'r') as f:
        chunkEnd = f.tell()
    while True:
        chunkStart = chunkEnd
        f.seek(size,1)
        f.readline()
        chunkEnd = f.tell()
        yield chunkStart, chunkEnd - chunkStart
        if chunkEnd > fileEnd:
            break

#init objects
pool = mp.Pool(cores)
jobs = []

#create jobs
for chunkStart,chunkSize in chunkify("input.txt"):
    jobs.append( pool.apply_async(process_wrapper,(chunkStart,chunkSize)) )

#wait for all jobs to finish
for job in jobs:
    job.get()

#clean up
pool.close()


# Another simple version ### This will read entire file into memory
import time
start_time = time.time()
from multiprocessing import Pool

def process_line(line):
    final_status =[]
    
    each_line = line.split('\t')
       
    if '@' not in each_line[0]:
        status = int(each_line[1])%16
        final_status.append(status)
    
    return final_status     
    #return "FOO: %s" % line
    #return line.split('\t')[1]



if __name__ == "__main__":
    pool = Pool(4)
    results=[]
    with open('/Volumes/Research_Data/MAC_Research_Data/Pan_cancer/CMT-SRR7780976-test2.sam') as source_file:
        # chunk the work into batches of 4 lines at a time
        results = pool.map(process_line, source_file)
print("--- %s seconds ---" % (time.time() - start_time))
print (results)

### 
import time
start_time = time.time()
def process_line(line):
    each_line = line.split('\t')
    print(each_line)
    """
    if '@' not in each_line[0]:
        status = int(each_line[1])%16
        final_status.append(status)
    
    return final_status
    """
if __name__ == "__main__":
    results =[]
    final_result=[]
    with open('/Volumes/Research_Data/MAC_Research_Data/Pan_cancer/CMT-SRR7780976-test2.sam') as source_file:
        for i in source_file :
            results = list(map(process_line, i))
            final_result.append(results)
            
            
            #results = list(map(process_line, i))
            #final_result.append(results)
        
print("--- %s seconds ---" % (time.time() - start_time))            

    


##### Read files line by lines

from multiprocessing import Process, Manager
import time
import itertools 

def do_work(in_queue, out_list):
    while True:
        item = in_queue.get()
        line_no, line = item

        # exit signal 
        if line == None:
            break

        # fake work
        time.sleep(.5)
        result = (line_no, line)

        out_list.append(result)


if __name__ == "__main__":
    num_workers = 4

    manager = Manager()
    results = manager.list()
    work = manager.Queue(num_workers)

    # start for workers    
    pool = []
    for i in range(num_workers):
        p = Process(target=do_work, args=(work, results))
        p.start()
        pool.append(p)

    # produce data
    with open("/Volumes/Research_Data/MAC_Research_Data/Pan_cancer/test.sam") as f:
        iters = itertools.chain(f, (None,)*num_workers)
        for num_and_line in enumerate(iters):
            work.put(num_and_line)

    for p in pool:
        p.join()

    # get the results
    # example:  [(1, "foo"), (10, "bar"), (0, "start")]
    print (sorted(results))


#read entire file into memory
import multiprocessing
from textwrap import dedent
from itertools import zip_longest

def process_chunk(d):
	"""Replace this with your own function
	that processes data one line at a
	time"""

	d = d.strip() + ' processed'
	return d 

def grouper(n, iterable, padvalue=None):
	"""grouper(3, 'abcdefg', 'x') -->
	('a','b','c'), ('d','e','f'), ('g','x','x')"""

	return zip_longest(*[iter(iterable)]*n, fillvalue=padvalue)

if __name__ == '__main__':

	# test data
	test_data = """\
	1 some test garbage
	2 some test garbage
	3 some test garbage
	4 some test garbage
	5 some test garbage
	6 some test garbage
	7 some test garbage
	8 some test garbage
	9 some test garbage
	10 some test garbage
	11 some test garbage
	12 some test garbage
	13 some test garbage
	14 some test garbage
	15 some test garbage
	16 some test garbage
	17 some test garbage
	18 some test garbage
	19 some test garbage
	20 some test garbage"""
	test_data = dedent(test_data)
	test_data = test_data.split("\n")

	# Create pool (p)
	p = multiprocessing.Pool(4)

	# Use 'grouper' to split test data into
	# groups you can process without using a
	# ton of RAM. You'll probably want to 
	# increase the chunk size considerably
	# to something like 1000 lines per core.

	# The idea is that you replace 'test_data'
	# with a file-handle
	# e.g., testdata = open(file.txt,'rU')

	# And, you'd write to a file instead of
	# printing to the stout

	for chunk in grouper(10, test_data):
		results = p.map(process_chunk, chunk)
		for r in results:
			print (r) 	# replace with outfile.write()
            