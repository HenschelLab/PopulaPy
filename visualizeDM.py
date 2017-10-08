
"""
## convert to binary
 2069  plink --file ALS16797 --make-bed --out ALS16797
 2071  plink --file alsafar2 --make-bed --out alsafar2
 2075  R
## identifying common snsp in R, for clean merger
https://www.biostars.org/p/56640/
 2080  plink  --bfile alsafar2 --extract list.snps --make-bed --out alsafar2c
 2081  plink  --bfile ALS16797 --extract list.snps --make-bed --out ALS16797c
## Merge
 2082  plink --bfile ALS16797c --bmerge alsafar2c.bed alsafar2c.bim alsafar2c.fam --out merged

## Working only with SNPs that have at least 5% allele frequency 
 2026  plink --file alsafar2 --maf 0.05 --make-bed --out alsafarMAF5
 2064  plink --file ALS16797 --maf 0.05 --make-bed --out ALS16797MAF5
 2085  plink --bfile ALS16797MAF5 --bmerge alsafarMAF5.bed alsafarMAF5.bim alsafarMAF5.fam --out mergedMAF5

## Creating Distance matrix 
 2086  plink --bfile mergedMAF5 --distance square --out mergedMAF5

## Created files: 
 1.1G  merged.bed
 149M  merged.bim
  19K  merged.fam
 6.9M  merged.log
 
 475M  mergedMAF5.bed
  65M  mergedMAF5.bim
  19K  mergedMAF5.fam
 1.2K  mergedMAF5.log
  10K  mergedMAF5.nosex

## Distance matrix  
  12M  mergedMAF5.dist
  10K  mergedMAF5.dist.id


MISC: creating homozygosity scores  
  plink --bfile mergedMAF5 --homozyg

## adding 3rd cohort file (only 1M SNPs as opposed to 5M)
plink --file CVRL --make-bed --out CVRL
  
"""

import pandas as pd
from scipy.spatial.distance import squareform
import scipy.cluster.hierarchy as sch
import pylab
import numpy as np
import pdb

from plinkio import plinkfile

dendrocolors = {}


def colordendrogram(Y, pos):
    def subcluster(subcl):
        subcl=int(subcl)
        if subcl < len(Y) + 1: ## singleton
            region = population(subcl).split('---')[-1].strip()
            dendrocolors[subcl] = popColorDict[region]
            return region
        else:
            return colordendrogram(Y, subcl - (len(Y)+1) )
    clusterId = pos + len(Y) + 1
    cl1, cl2, tm, clusterSize =  Y[pos]
    result1 = subcluster(cl1)
    result2 = subcluster(cl2)
    if (result1 == 'Mixed' or result2 == 'Mixed') or (result1 != result2):
        dendrocolors[clusterId] = 'gray'
        return 'Mixed'
    else:
        dendrocolors[clusterId] = popColorDict[result1]
        return result1

popDict = {('America', 'America'): 'America',
           ('Asia', 'Central_South_Asia'): 'CentralSouthAsia',
           ('Asia', 'Est_Asia'): 'EastAsia',
           ('Europe', 'Europe'): 'Europe',
           ('Middle_Est', 'Middle_Est'): 'MiddleEast',
           ('North_Africa', 'Middle_Est'): 'NorthAfrica',
           ('Oceania', 'Oceania'): 'Oceania',
           ('Subsaharian_Africa', 'Africa'): 'Subsahara'}

popColorDict = {'Oceania': 'g',
                'Europe': 'y',
                'America': 'brown',
                'EastAsia': 'indigo',
                'CentralSouthAsia': 'm',
                'UAE':'lightgreen',
                'MiddleEast': 'c',
                'NorthAfrica': 'b',
                'Subsahara': 'r',
                'Mixed': 'gray' }

def population(nid):
    sampleId = ids.iloc[nid] ## sampleId 0 -> HGDP00448 -> 
    if sampleId in sampleInfo.index:
        row = sampleInfo.loc[sampleId]
        return "%s/%s/%s   ---  %s"%(row['population'], row['Pop7Groups'],  row['Sex'], popDict[(row['Region'], row['Pop7Groups'])])
    else:
        return str(sampleId) + " --- UAE"


def populationSelectedReps(nid):
    sampleId = ids.iloc[nid] ## sampleId 0 -> HGDP00448 ->
    if sampleId in sampleInfo.index: ## HGDP dude
        row = sampleInfo.loc[sampleId]
        return "" #popDict[(row['Region'], row['Pop7Groups'])][:2]
    else:
        if str(sampleId) in reps:
            return str(sampleId)
        else:
            #if not sampleId.startswith("HGDP"):
            #    pdb.set_trace()
            return ""

    
droppingRelatives = False

basename = 'uae_hgdp1LD'
    
wdir = "/research/gutsybugs/KUMI/Data/"
hetfile = "%s.het" % basename ## to be calculated
plink_file = plinkfile.open(basename)

dmfile = "%s.mdist" % basename  #"mergedQCLD.mdist" # "merged_1ibs.mdist" #"merged.dist"
dmfileIDs = dmfile + ".id"
## Loading distance matrix, produced by plink, putting ids to dataframe
ids = pd.read_csv(dmfileIDs, delimiter='\t', header=None)[1]
dm = pd.read_csv(dmfile, header=None, delimiter='\t')
het = pd.read_csv(hetfile, delimiter='\s+')
sampleInfo = pd.read_csv("hgdp/HGDPid_populations.csv", sep=',',  index_col='Id')
fam = pd.read_csv('%s.fam' % basename, sep=' ', header=None)

#reps = map(str, set(pd.read_csv('100medoids.txt', header=None)[0]))
reps = '10187 12742 13120 13076 10651 10347 10215 10725 10926 12599'.split()

#selectSamples = 100
if basename == 'uae_hgdp1LD':
    selectIDs = pd.read_pickle("%s_selection.pcl"% basename)
    dm = dm.iloc[selectIDs, selectIDs]
    ids = ids[selectIDs]
    dm.index = dm.columns = ids
else:
    dm.index = dm.columns = ids

sample_list = plink_file.get_samples( )
locus_list = plink_file.get_loci( )

#makeSampleID = lambda sample: "%s_%s"%(sample.fid, sample.iid)
sampleIDs = [sample.iid for sample in sample_list]

makerow = lambda row,locus: [snp for snp in row] + [locus.chromosome, locus.name, locus.position, locus.bp_position]

#df = pd.DataFrame(data=(makerow(row, locus) for locus, row in zip(locus_list, plink_file)),
#                  index=np.arange(len(locus_list)),
#                  columns=sampleIDs + ['chr', 'name', 'genDist', 'pos'])

print "Loading done"

relThreshold = 0.08
## removing relatives, based on distance histogram, see results
if droppingRelatives:
    x,y = np.where(dm.values < relThreshold)
    dropindex = ids[[a for a,b in zip(x,y) if a!=b]]
    print "Dropping %s (related) people" %len(dropindex)
    dm = dm.drop(dropindex,axis=1).drop(dropindex)

## scale the rest, like Z-score but only positive
scaleDists = False
if scaleDists:
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

fya = 0.01

fig = pylab.figure(figsize=(35., 10))

## DENDROGRAM
print "calculating dendrogram for HC"
axdendro = fig.add_axes([fx1+mrg,0.6,fx2-fx1-mrg,0.3])
Y = sch.linkage(squareform(dm), method='complete') ## consider NJ!
colordendrogram(Y, len(Y)-1)
Z = sch.dendrogram(Y, leaf_label_func=populationSelectedReps, link_color_func = lambda k: dendrocolors[k])
print "done"

index = dm.columns[Z['leaves']]
dm = dm[index].reindex(index) ## reorder dm according to dendrogram
cmap = pylab.cm.jet

## DISTANCE MATRIX

if False:
    axmatrix = fig.add_axes([fx1+mrg,mrg, fx2-fx1-mrg, 0.5])
    im = axmatrix.matshow(dm.values, aspect='auto', cmap=cmap)

## HETEROZYGOSITY SCORE
#cmap2 = pylab.cm.bwr
#axmatrix2 = fig.add_axes([fx1+mrg, 0.01, fx2-fx1-mrg, 0.02])
#im2 = axmatrix2.matshow( [het.iloc[index].F], aspect='auto', cmap=cmap2)

## ADMIXTURE BARPLOT
## expects the plot for 8 populations (ideally using supervised)
if False:
    k = 8
    q = pd.read_csv('%s.%s.Q'% (basename, k), header=None, sep=' ')
    q['Id'] = fam[1]
    q.set_index('Id', inplace=True)
    q1 = q.loc[index]
    axmatrix3 = fig.add_axes([fx1+mrg, 0.25, fx2-fx1-mrg, 0.3])
    axmatrix3.set_ylim([0,1])
    colors = ['r', 'b', 'c', 'm', '#624ea7', 'g', 'yellow', 'k', 'maroon']
    q1.plot.bar(stacked=True, linewidth=0, width=1, ax=axmatrix3, color = colors[:k], legend=False)

#q1 = q[index]
#cmap3 = pylab.cm.binary
#axmatrix3 = fig.add_axes([fx1+mrg, 0.6+mrg, fx2-fx1-mrg, 0.15])
#df1 = df[df.chr==1].iloc[:1500]
#im3 = axmatrix3.bar(q, aspect='auto', cmap=cmap3)


#cmap3 = pylab.cm.binary
#axmatrix3 = fig.add_axes([fx1+mrg, 0.02+mrg, fx2-fx1-mrg, fy1-(0.02+mrg)])
#df1 = df[df.chr==1].iloc[:1500]
#im3 = axmatrix3.matshow( df1[index], aspect='auto', cmap=cmap3)


#fig.show()
#fig.savefig('%s_admix3a.svg' % basename)
#fig.savefig('%s_admix3a.png' % basename)
