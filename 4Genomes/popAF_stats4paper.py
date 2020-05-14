import pandas as pd
import numpy as np
from collections import Counter

"""                                                                                                                                                                                                                                                                            
allele frequencies (and their loci) annotated with SnpEff, clinvar, ccdsGene and Kegg using vtools,                                                                                                                                                                            
applied to the table of                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                               
generates numbers for table 1 in Population paper                                                                                                                                                                                                                              
"""

if __name__ == "__main__":
    afSheet = "../Data/z4_eff_cln.csv.gz" ##                                                                                                                                                                                                                                   
    af = pd.read_csv(afSheet, sep='\t')
    ## embellishing the table a bit                                                                                                                                                                                                                                            
    af.columns = [col.strip() for col in af.columns.values]
    for col in af.columns.values: ## making proper floats (possibly doable during import/dtype specs)                                                                                                                                                                          
        if col.startswith('AF_'):
            af[col] = af[col].apply(lambda item: np.nan if type(item)==str and item.strip()=='.' else float(item))

    for col in 'EFF CLNSIG KgDesc ccdsGene_name'.split():
        af[col] = af[col].str.strip()
    ## SnpEFF comes like this: INTRON(MODIFIER||||184|LCN15|protein_coding|CO..., first part (INTRON) to new col EFF1                                                                                                                                                          
    af['EFF1'] = af.apply(lambda row: row['EFF'].split('(')[0], axis=1)


d1 = Counter(af[af['AF_Z']> 4]['EFF1'])
d2 = Counter(af[af['AF_Z']<-4]['EFF1'])
df = pd.DataFrame([d2, d1])

df['Others'] = df.UPSTREAM + df.DOWNSTREAM + df.SYNONYMOUS_CODING + df.UTR_3_PRIME + df.UTR_5_PRIME + df.INTERGENIC

cols = ['NON_SYNONYMOUS_CODING', 'START_GAINED',
        'STOP_GAINED',  'EXON', 'INTRON', 'SPLICE_SITE_REGION', 'Others']
df1 = df[cols]
df1['Total'] = df1.sum(axis=1)
print(df1.append(df1.sum(numeric_only=True), ignore_index=True))