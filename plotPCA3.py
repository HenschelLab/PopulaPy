"""
Mixing Several ScatterPie plots

Example script for using scatter pie.
Draws three types of scatter pies:
1. selection 0: selected representatives are highlighted - red edge
2. selection 1: remainders (non-representatives from UAE study) - green edge
3. selection 2: Context: HGDP samples - smaller circles, no edge, 80% opac

Requires * eigenvec file to be generated with plink
         * fam file from running plink of combined dataset 
         * Outputfile from admixture (.Q)         
"""
import matplotlib
matplotlib.use('agg')

import pandas as pd
import matplotlib.pyplot as plt
from scatterPie import scatterpie

# plink -bfile uae_hgdp1LD --keep sel5.txt --make-bed --out uae_hgdp5pop
hgdpSampleInfo = pd.read_csv("hgdp/HGDPid_populations.csv", sep=',', index_col='Id')
k = 8
## Taking admix data from full run
## paste uae_hgdp1LD.8.Q uae_hgdp1LD.fam > uae_hgdp1LDa.8.Q
q = pd.read_csv('uae_hgdp1LDa.8.Q', header=None, sep='\s+', index_col=9)
qcol = ['Com%s'%i for i in range(k)]
q.columns = qcol + ['Fam', 'P1', 'P2', 'P3', 'P4']
colors = ['r','b','c','m','indigo','g','y', 'saddlebrown']
    
for ri in range(10):
    basename = f'uae_hgdp5pop_125_{ri}'
    fam = pd.read_csv('%s.fam' % basename, sep=' ', header=None, index_col=1) 
    ev = pd.read_csv('%s.eigenvec'%basename, header=None, sep='\s+', index_col=1) ## assuming this is arranged in the correct order
    ev.columns = ['Fam2'] + ['PC%02d'%i for i in range(1,21)]

    j = pd.concat([fam, ev, q], join='inner', axis=1)
    hgdp = j[j.index.str.startswith('HGDP')]
    uae = j[~j.index.str.startswith('HGDP')]


    fig = plt.figure()
    ax = fig.add_subplot(111) 

    scatterpie(hgdp[['PC01','PC02']].values, hgdp[qcol].values, [3]*len(hgdp), colors=colors, alpha=0.8, ax=ax, lw=0) ## HGDP
    scatterpie(uae[['PC01','PC02']].values, uae[qcol].values, [10]*len(uae), colors=colors, ax=ax, edgecolor='black', lw=0.1) #, labels=list(fam[selection][1]))

    plt.savefig('../Results/%s_scatterPie.svg' %basename)
    plt.savefig('../Results/%s_scatterPie.png' %basename)
    plt.savefig('../Results/%s_scatterPie.pdf' %basename)
    print("finished", basename)
    plt.clf()

