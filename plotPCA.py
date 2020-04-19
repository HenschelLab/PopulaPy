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

popColorDict = {'Oceania': 'g',
                'Europe': 'y',
                'America': 'k',
                'EastAsia': 'crimson',
                'CentralSouthAsia': 'm',
                'UAE':'lightgreen',
                'MiddleEast': 'c',
                'NorthAfrica': 'b',
                'Subsahara': 'r',
                'Palestinian': 'darkseagreen',
                'Druze': 'Lime',
                'Beduin': 'Cyan',
                'Mixed': 'gray' }
#popDict = {('America', 'America'): 'America',
#           ('Asia', 'Central_South_Asia'): 'CentralSouthAsia',
#           ('Asia', 'Est_Asia'): 'EastAsia',
#           ('Europe', 'Europe'): 'Europe',
#           ('Middle_Est', 'Middle_Est'): 'MiddleEast',
#           ('North_Africa', 'Middle_Est'): 'NorthAfrica',
#           ('Oceania', 'Oceania'): 'Oceania',
#           ('Subsaharian_Africa', 'Africa'): 'Subsahara'}

#def assignColor(sampleId):
#    if sampleId in hgdpSampleInfo.index:
#        row = hgdpSampleInfo.loc[sampleId]
#        return popColorDict[popDict[(row['Region'], row['Pop7Groups'])]]
#    else:
#        return popColorDict['UAE']

basename = 'mergedQC_CVRL_flip_hgdp2' ## combined file 'uae_hgdp1LD'
hgdpSampleInfo = pd.read_csv("hgdp/HGDPid_populations.csv", sep=',', index_col='Id')
fam = pd.read_csv('%s.fam' % basename, sep=' ', header=None) 
ev = pd.read_csv('%s.eigenvec'%basename, header=None, sep='\s+') ## assuming this is arranged in the correct order

#colors = [assignColor(sid) for sid in ev[1]]
k = 8
q = pd.read_csv('%s.%s.Q'% (basename, k), header=None, sep=' ')
#q['Id'] = fam[1]
#q.set_index('Id', inplace=True)
#reps = '13050 10898 12906 12982 12677 10596 10288 11045 10177 10416 10231 11241 120383 10300 13195'.split()
#reps1 = map(str, range(1,139))
#reps2 = map(str, range(139,202))
#reps3 = map(str, range(202,274))
## Mark selected representatives
reps = '12584 10745 13024 120075'.split()
#reps = ['10885', '12553', '12657', '12960', '12964', '13024', '13049', '13116', '13119', '120060', '120067', '120069', '120075', '120129', '120131', '120133', '120136', '120164', '120165', '120190', '120195', '120197', '120213', '120287', '12608', '120233'] # outsourced seq.
#reps = ['12742', '9239299151_R03C01']
#reps = '10187 12742 13120 13076 10651 10347 10215 10725 10926 12599'.split() ## those sent for sequencing
#reps = '10187 10215 10725 10926 12599'.split()
#reps = map(str, list(pd.read_csv('100selected.csv')['Sample ID']))
## seriously ugly code, in need of a rewrite
selection0 = [f in reps for f in fam[1]] ## True/False boolean mask
#selection_kw1 = [f in reps1 for f in fam[1]] ## True/False boolean mask
#selection_kw2 = [f in reps2 for f in fam[1]] ## True/False boolean mask
#selection_kw3 = [f in reps3 for f in fam[1]] ## True/False boolean mask
selection1 = [not f.lower().startswith('hgdp') and not f in reps for f in fam[1]] ## True/False boolean mask #and not f in reps
selection2 = [f.lower().startswith('hgdp') for f in fam[1]] ## True/False boolean mask

#ev_kw1 = ev[selection_kw1]
#q_kw1 = q[selection_kw1]
#ev_kw2 = ev[selection_kw2]
#q_kw2 = q[selection_kw2]
#ev_kw3 = ev[selection_kw3]
#q_kw3 = q[selection_kw3]
ev0 = ev[selection0]
q0 = q[selection0]
ev1 = ev[selection1]
q1 = q[selection1]
ev2 = ev[selection2]
q2 = q[selection2]

fig = plt.figure()
ax = fig.add_subplot(111) 
#colors10 = ['r','b','c','darkseagreen','Lime','m','indigo','g','y', 'saddlebrown'] ## use for 10 ancestral pops (*.10.Q)
colors = ['r','b','c','m','indigo','g','y', 'saddlebrown']

scatterpie(zip(ev2[2], ev2[3]), q2.values[:, :k-1], [5]*len(q2), colors=colors, alpha=0.8, ax=ax, lw=0) ## HGDP
scatterpie(zip(ev1[2], ev1[3]), q1.values[:, :k-1], [10]*len(q1), colors=colors, ax=ax, edgecolor='black', lw=0.15) #, labels=list(fam[selection][1])) #UAE, non-representatives
#scatterpie(zip(ev1[2], ev1[3]), q1.values[:, :k-1], [3]*len(q1), colors='c', ax=ax, lw=0.00) # version for 1 genome paper
#representatives
scatterpie(zip(ev0[2], ev0[3]), q0.values[:, :k-1], [30]*len(q0), colors=colors, ax=ax, edgecolor='lightgreen', lw=0.15) # labels=list(fam[selection0][1])) 
# 100 reps#scatterpie(zip(ev_kw1[2], ev_kw1[3]), q_kw1.values[:, :k-1], [5]*len(q_kw1), colors=colors, ax=ax, edgecolor='red', lw=0.2)  ## Kuwait 1
#scatterpie(zip(ev_kw2[2], ev_kw2[3]), q_kw2.values[:, :k-1], [5]*len(q_kw2), colors=colors, ax=ax, edgecolor='black', lw=0.2) ## Kuwait 2
#scatterpie(zip(ev_kw3[2], ev_kw3[3]), q_kw3.values[:, :k-1], [5]*len(q_kw3), colors=colors, ax=ax, edgecolor='green', lw=0.2) ## Kuwait 2

basename = '20outsourcedGenomes1' 
plt.savefig('../Results/%s_scatterPie.svg' %basename)
plt.savefig('../Results/%s_scatterPie.png' %basename)
plt.savefig('../Results/%s_scatterPie.pdf' %basename)

#plt.show()

#plt.scatter(ev[2], ev[3], color=colors)
