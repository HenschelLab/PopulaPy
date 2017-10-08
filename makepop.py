import pandas as pd
import matplotlib.pyplot as plt
## Define 8 populations:

popDict = {('America', 'America'): 'America',
           ('Asia', 'Central_South_Asia'): 'CentralSouthAsia',
           ('Asia', 'Est_Asia'): 'EastAsia',
           ('Europe', 'Europe'): 'Europe',
           ('Middle_Est', 'Middle_Est'): 'MiddleEast',
           ('North_Africa', 'Middle_Est'): 'NorthAfrica',
           ('Oceania', 'Oceania'): 'Oceania',
           ('Subsaharian_Africa', 'Africa'): 'Subsahara'}

def population(row):
    return popDict[(row['Region'], row['Pop7Groups'])]

uaeSNPs = "uae_hgdp1"
fam = pd.read_csv('%s.fam' % uaeSNPs, sep=' ', header=None)
    
ids = pd.read_csv("hgdp/HGDPid_populations.csv", sep=',')
ids.set_index('Id', inplace=True)

for id, row in fam.iterrows():
    if row[1] in ids.index:
        print population(ids.loc[row[1]])
    else:
        print '--'

        
    
