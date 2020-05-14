import pandas as pd
import numpy as np
import sys
import gzip

"""
Workflow:

1. parseGnomad: 
Input: allele frequency tables, Gnomad-vcf files (currently r2.0.2, better r2.1 -contains SAS, south asian)
Output: Data/popAF*csv allele frequency spreadsheat, combining Gnomad, GME*, UAE
(gzip popAF*csv)

2. frqAnalysis
popAFs4_*.csv - calculates Z-scores, adds it to table

3. frqAnalysis_report, oo version of 2.,
also generates BED file for javascript based visualization (ideogram.js)

4. table2vcf
PROBLEM: Memory error for large chromosomes 1, 2. 3,6,7,8 got killed as well, possibly for the same reason

takes a pandas generated (,-sep) csv, creates a vcf
NEW: to save memory

"""

cmd = "nohup python table2vcf.py %s  > t2vcf.nohup%s.out  2> t2vcf.nohup%s.err &"
info = '''##INFO=<ID=AF_%s,Number=A,Type=Float,Description="Allele Freq for %s">'''
infoZ ='''##INFO=<ID=AF_gMAX,Number=A,Type=Float,Description="Max Gnomad Allele Freq Z-score">
##INFO=<ID=AF_gMIN,Number=A,Type=Float,Description="Min Gnomad Allele Freq Z-score">
##INFO=<ID=AF_Z,Number=A,Type=Float,Description="UAE Allele Freq Z-score">
##INFO=<ID=AF_Z1,Number=A,Type=Float,Description="UAE Allele Freq Z-score, also considering GME AF">
'''
pops = ['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH', 'GME', 'UAE']

meta = '''##fileformat=VCFv4.0
##fileDate=20190131
##reference=GRCh37.p13
'''

afs = 'AFR,AMR,ASJ,EAS,FIN,NFE,OTH,GME,UAE,gMAX,gMIN,Z,Z1'.split(',')
header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
if len(sys.argv) == 1:
    for i in [1,2,3,6,7,8]: #list(range(1,23)) + [25]:
        print (cmd % (i,i,i))
else:
    chrom = int(sys.argv[-1])
    with open("Data/popAFs4_%s.vcf"%chrom, "w") as w:
        ## header
        w.write(meta)
        for pop in pops:
            w.write(info%(pop,pop)+"\n")
        w.write(infoZ)
        w.write(header)
        df = pd.read_csv('Data/popZ4AFsFull_%s.csv.gz'%chrom)
        for line in gzip.open('Data/popZ4AFsFull_%s.csv.gz'%chrom):
            if line.decode().startswith(',CHROM'): continue
            row=line.decode().strip().split(',')
            bla=";".join(["AF_%s=%s"% (afs[i],v) for i,v in enumerate(row[8:]) if v])
            meta=[row[i] for i in [1,2,5,3,4]] 
            w.write('\t'.join(meta + ['.', '.'] + [bla]) + "\n")
    print('Done with chrom. %s' % chrom)