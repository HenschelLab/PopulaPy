import pandas as pd
import subprocess
import numpy as np

hgdpSampleInfo = pd.read_csv("hgdp/HGDPid_populations.csv", sep=',', index_col='Id')
popDict = {('America', 'America'): 'America',
           ('Asia', 'Central_South_Asia'): 'CentralSouthAsia',
           ('Asia', 'Est_Asia'): 'EastAsia',
           ('Europe', 'Europe'): 'Europe',
           ('Middle_Est', 'Middle_Est'): 'MiddleEast',
           ('North_Africa', 'Middle_Est'): 'NorthAfrica',
           ('Oceania', 'Oceania'): 'Oceania',
           ('Subsaharian_Africa', 'Africa'): 'Subsahara'}
from collections import Counter

f = lambda row: popDict[(row['Region'], row['Pop7Groups'])]
hgdpSampleInfo['ContPop'] = hgdpSampleInfo.apply(f, axis=1)

basename = 'uae_hgdp1LD'
fam = pd.read_csv('%s.fam' % basename, sep=' ', header=None, index_col=1)
fam.index.name = 'Id'
#famH = pd.concat([fam, hgdpSampleInfo], axis=1, sort=True)
famH = pd.merge(fam, hgdpSampleInfo, how='left', on='Id')
famH['ContPop'].fillna('UAE', inplace=True)

sampleStat = {'UAE': 1000,
              'CentralSouthAsia': 210,
              'Subsahara': 127,
              'Oceania': 39,
              'Europe': 161,
              'MiddleEast': 148,
              'America': 108,
              'EastAsia': 241,
              'NorthAfrica': 30}

## select only relevant populations (major contributoris)
#sel = famH[famH['ContPop'].isin(['CentralSouthAsia', 'UAE', 'Subsahara', 'MiddleEast', 'NorthAfrica', 'Europe'])]
## subsampling
selSize = 125

for ri in range(10):
    selection = []
    for pop, size in sampleStat.items():
        if pop in 'Oceania America EastAsia NorthAfrica'.split(): continue
        s = list(np.random.choice(famH[famH.ContPop==pop].index.values, selSize))
        selection += s

    famH2 = famH.loc[selection]
    randSelFile = 'sel5_%s_%s.txt' % (selSize,ri)
    w = open(randSelFile, 'w')
    for i, row in famH2.iterrows():
        if np.isnan(row[0]) or row.ContPop in 'Oceania America EastAsia NorthAfrica'.split(): continue
        print(int(row[0]), i, file = w)
    w.close()
    cmd = f"plink -bfile uae_hgdp1LD --keep sel5_{selSize}_{ri}.txt --make-bed --out uae_hgdp5pop_{selSize}_{ri}"
    cmd2 = f"plink -bfile uae_hgdp5pop_{selSize}_{ri} --pca --out uae_hgdp5pop_{selSize}_{ri}"
    print(cmd)
    print(cmd2)
    #subprocess.call(cmd, shell=True)
