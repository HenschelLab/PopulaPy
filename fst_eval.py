import glob
import pandas as pd
import subprocess
import numpy as np

def avgFst(i,j):
    fst = pd.read_csv(f'FstUAE_HGDP/pairFst_{i}_{j}.weir.fst', sep='\t') 
    return fst.WEIR_AND_COCKERHAM_FST.mean()

def procline(line):
    job = int(line.split('.')[1]) -1
    val = float(line.strip().split()[-1])
    return job, val
    
pops = sorted(glob.glob('uaeIDs_*.txt')) + sorted(glob.glob('hgdp_*.txt'))
fst = np.zeros((len(pops), len(pops)))

## Reading pairwise
#fstLog = dict(procline(line) for line in open('weightedFst04.log'))
#fstLogDiagonal = dict(procline(line) for line in open('weightedFst99.log'))
fstLog = dict(procline(line) for line in open('FstUAE_HGDP/meanFst04.log'))
fstLogDiagonal = dict(procline(line) for line in open('FstUAE_HGDP/meanFst99.log'))

    
for idx, (i,j) in enumerate([(k,l) for k in range(len(pops)) for l in range(k)]):
    fst[i,j] = fst[j,i] = fstLog[idx] #avgFst(i,j)

for i in range(len(pops)):
    fst[i,i] = fstLogDiagonal[i]
    print (pops[i], fst[i,i])

prettyPops = [pop[0].upper() + '_'.join(pop.split('_')[1:])[:-4] for pop in pops]
F = pd.DataFrame(fst, index=prettyPops, columns=prettyPops)
F.to_pickle('pairwiseMeanFst.pcl')
