#!/usr/bin/python3

### python DbSNP_vcf_filtering.py *.vcf <dbSNP.vcf>

import sys
import glob

def dbsnp_dic(dbsnp):
    result = {}
    with open(dbsnp,'r') as f:
        for lst in f:
            if lst[0] != '#':
                info = lst.split('\t')
                chrom = info[0]
                pos = info[1]
                ref = info[3]
                alt = info[4]
                Key = chrom+":"+pos
                value = ref+":"+alt
                if Key not in result.keys():
                    result[Key] = [value]
                else:
                    result[Key].append(value)

    return result


file_name = sys.argv[1]
dbsnp = sys.argv[2]
# dbsnp dic
dbsnpDic = dbsnp_dic(dbsnp)
# filter vcf files
with open(file_name,'r') as f:
    file = f.read()

out = open(file_name.split('.vcf')[0] + 'Line_to_line_tuple_dbSNPfiltered.vcf', 'w')
lst = file.split('\n')[:-1]
for i in range(len(lst)):
    if lst[i][0] == '#':
        out.write(lst[i] + '\n')
    else:
        info = lst[i].split('\t')
        chrom = info[0]
        pos = info[1]
        ref = info[3]
        alt = info[4]
        Key = chrom+":"+pos
        Value = ref+":"+alt
        #status = info[6]
        if Key not in dbsnpDic.keys():
            out.write(lst[i] + '\n')
        else:
            if Value not in dbsnpDic[Key]:
                out.write(lst[i] + '\n')

       
out.close()

