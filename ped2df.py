from plinkio import plinkfile
import pandas as pd
import numpy as np
import MySQLdb
from MySQLdb.cursors import DictCursor
from scipy.spatial.distance import squareform
import scipy.cluster.hierarchy as sch
import pylab

#filename = '/research/gutsybugs/Software/Plink/Tutorial/hapmap1'
wdir = "/research/gutsybugs/KUMI/Data/"
filename = '%s/mergedQC' % wdir
dmfile = "%s.dist" % filename
hetfile = "%s.het" % filename
plink_file = plinkfile.open(filename)
#plink_file = plinkfile.open('/research/gutsybugs/KUMI/Data/mergedQC') ## TAKES LOOOOOONG!

sample_list = plink_file.get_samples( )
locus_list = plink_file.get_loci( )

makeSampleID = lambda sample: "%s_%s"%(sample.fid, sample.iid)
sampleIDs = [makeSampleID(sample) for sample in sample_list]

ntdict = {'N': np.NaN, 'A': 1, 'C': 2, 'G': 3, 'T': 4}
if filename.endswith('hapmap1') or filename.endswith('mergedQC'):
    snpFct = lambda x: x
else:
    snpFct = lambda x: ntdict[x]
    
makerow = lambda row,locus: [snpFct(snp) for snp in row] + [locus.chromosome, locus.name, locus.position, locus.bp_position]

df = pd.DataFrame(data=(makerow(row, locus) for locus, row in zip(locus_list, plink_file)),
                  index=np.arange(len(locus_list)),
                  columns=sampleIDs + ['chr', 'name', 'genDist', 'pos'])


dmfileIDs = dmfile + ".id"
## Loading distance matrix, produced by plink, putting ids to dataframe
ids = pd.read_csv(dmfileIDs, delimiter='\t', header=None)[0] ## !!! Family ID (for HapMap only)
dm = pd.read_csv(dmfile, header=None, delimiter='\t')
het = pd.read_csv(hetfile, delimiter='\s+') ## Those idiots!! Dont switch between csv dialects, be consistent!

relThreshold = 0.08

dms = squareform(dm)
dms = dms[dms > relThreshold]
minv = (dms.min()-dms.mean())/dms.std()
scalefct = lambda x: max(0, (x-dms.mean())/dms.std() - minv)
#scalefct = lambda x: max(0, x - dms.min())
dm = dm.applymap(scalefct)
print "scaled distances"

fx1 = 0.05
fx2 = 0.95
mrg = 0.01
fy1 = .7
fy2 = .85

fya = 0.01

fig = pylab.figure(figsize=(25., 40))
axdendro = fig.add_axes([fx1+mrg,fy2+.03,fx2-fx1-mrg,0.11])

Y = sch.linkage(squareform(dm), method='complete') ## consider NJ!
Z = sch.dendrogram(Y)

index = dm.columns[Z['leaves']]
dm = dm[index].reindex(index) ## reorder dm according to dendrogram
cmap = pylab.cm.jet

axmatrix = fig.add_axes([fx1+mrg,fy1+mrg, fx2-fx1-mrg, fy2-(fy1+mrg)])
im = axmatrix.matshow(dm.values, aspect='auto', cmap=cmap)

cmap2 = pylab.cm.bwr
axmatrix2 = fig.add_axes([fx1+mrg, 0.01, fx2-fx1-mrg, 0.02])
im2 = axmatrix2.matshow( [het.iloc[index].F], aspect='auto', cmap=cmap2)

cmap3 = pylab.cm.binary
axmatrix3 = fig.add_axes([fx1+mrg, 0.02+mrg, fx2-fx1-mrg, fy1-(0.02+mrg)])
im3 = axmatrix3.matshow( df[df.chr==1].iloc[:1500, index], aspect='auto', cmap=cmap3)

fig.savefig("hapmap.png")
