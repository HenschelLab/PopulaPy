import pandas as pd
import matplotlib.pyplot as plt

uaeSNPs = "uae_hgdp1"
populations = 7

ids = pd.read_csv("HGDPid_populations.csv", sep=',')
ids.set_index('Id', inplace=True)
#order = [line.split()[1] for line in open("uae_hgdp1.fam")]
fam = pd.read_csv('hgdp1.fam', sep=' ', header=None)

## plotting ADMIXTURE barplot
q = pd.read_csv('hgdp1.7.Q', header=None, sep=' ')
q['Id'] = fam[1]
q.set_index('Id', inplace=True)
ids1 = ids.sort_values(by=['Pop7Groups', 'Region', 'Geographic_origin', 'population', 'Sex'])
q1 = q.loc[ids1.index.values] ## sort admixture barplot by population etc

#q1.plot.bar(stacked=True, linewidth=0, width=0.95)

#plt.show()

## hgdp = pd.read_csv("HGDP_FinalReport_Forward.txt", header=None, sep='\t')
## sexDict = {'M': '1', 'F': '2'}
## hgdp.set_index(0)

## w = open("hgdp.ped", 'w')
## for colID in hgdp.columns.values[1:]:
##     col = hgdp.loc[:, colID]
##     sampleID = col.loc[0]
##     sample = ids.loc[sampleID]
##     sex = sexDict.get(sample['Sex'], '0')
##     print >> w, "\t".join([str(colID), sampleID, '0', '0', sex, '0']) + '\t' + "\t".join(["%s\t%s" % (geno[0], geno[1]) for geno in col.loc[1:]])
## w.close()
    
