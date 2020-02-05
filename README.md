# QC_for_NGS_data

The target_reads.py script can analyze how many sequence reads are in the exon regions that were definfed in the exon interval files to see the whole exome sequence quality. 

Idealy, the sequence file needs to have at least 70% reads in the exon regions, but if we use CDS region, then it better has > 50% reads

The script take sequence sam file , exon interval file and sequence reads as input files.

The script implements the binary search to reduce the searching time. 

One whole exome sequence file takes about 1.5 hours to get the exon(CDS) mapping rates using python script.
One WES file take about 15 mins to get the exon (CDS) mapping rates using java script