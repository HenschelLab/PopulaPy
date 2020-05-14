"""comparison of our allele frequencies with GnomAD and GME
## Algorithm:
Parse GnomAD's genome uncompressed vcf files (r2.0.2)
into 

Note: GnomAD underwent a pretty dramatic change from r2.0 -> r2.1
for parsing the INFO column for AFs, see annotate.py
"""
import pdb
import sys
import numpy as np
import pandas as pd
import re
from collections import defaultdict
from Bio import bgzf

cmd = "nohup python parseGnomad.py %s > Nohup/parseGnomad.%s.nohup.out 2> Nohup/parseGnomad.%s.nohup.err &"
popAF = re.compile("AF_[A-Z]{3}=")
wdir = "/bmshare/ahenschel/KUMI/Data"
gnomadDir = "/bmshare/gihan/GnomAD"
#gnomadDir = "/bmshare/ahenschel/References/GnomadGenomes"
gmeDir = "/bmshare/ahenschel/References/GME"

class VCFLine:
    def __init__(self, line):
        self.chrom, pos, self.rsid, self.ref, self.alt, _, _, self.info = line.strip().split()
        self.alts = self.alt.split(',')
        self.indel = len(self.alts[0]) > 1
        self.pos = int(pos)

class GnomADLine(VCFLine):
    def info2popfreqs(self, whichAlt):
        def getCorrectAllele(keyVals):
            key, vals = keyVals.split('=')
            val = vals.split(',')[whichAlt]
            try:
                valnum = float(val)
            except ValueError:
                valnum = np.nan
            return key, valnum 
        self.popfreqs = [getCorrectAllele(keyVals) for keyVals in self.info.split(';') if popAF.match(keyVals) ]

    def createPopfreqDF(self, gmeAF=None): # gmeFreqs['GME_AF']
        if not gmeAF is None:
            self.popfreqs.append(('AF_GME', gmeAF))
        self.af = pd.DataFrame.from_records(self.popfreqs)

    def closestMatch(self, maf):
        self.af['diff'] = abs(self.af[1] - maf)
        return self.af.sort_values('diff').iloc[0]

        
class Annotation: ## somewhat compatible with VCFLine
    def __init__(self, frqOmni):
        self.ref = frqOmni.ref
        self.alts = frqOmni.split(',') ## should be only one though (?)

class FreqStats:
    def __init__(self, frqOmni):
        if type(frqOmni) == pd.DataFrame:
            frqOmni = frqOmni.iloc[0] ## just taking the first, could also take one with higher NCHROBS, rsid over kgp etc...
        self.minorAllele = frqOmni.A1 ## A1: plink's minor allele (usually?!?!?!) Not reliable! odes not correspond to ref/alt in hg19 :,'-(
        self.majorAllele = frqOmni.A2
        self.pos = frqOmni._name
        self.maf = frqOmni.MAF

    def compatibilityCheck(self, annotation, db='GnomAD'):
        global correctMatch, flippedMatch
        """checking compatibility between GnomAD/GME variant and our omni"""
        if self.majorAllele == annotation.ref and self.minorAllele in annotation.alts: 
            self.whichAlt = annotation.alts.index(self.minorAllele)
            correctMatch += 1
            return True
        elif self.minorAllele == annotation.ref and self.majorAllele in annotation.alts: 
            tmp = self.minorAllele
            self.minorAllele = self.majorAllele
            self.majorAllele = tmp
            self.whichAlt = annotation.alts.index(self.minorAllele)
            self.maf = 1-self.maf
            flippedMatch += 1
            return True
        return False

class Annotation:
    def __init__(self, row):
        self.alts = [row.alt]
        self.ref = row.ref
        self.frq = row.GME_AF
        
def makeRow(gnomadLine, fs):
    uaeMinor = fs.minorAllele if fs else 'X'
    meta = [gnomadLine.chrom, gnomadLine.pos, gnomadLine.ref, gnomadLine.alt, gnomadLine.rsid, len(gnomadLine.alts), uaeMinor]
    freqs = [dict(gnomadLine.popfreqs).get('AF_%s'%pop, '') for pop in pops]    
    return sep.join(map(str, meta + freqs)) +'\n'

def lookupRS(iid):
    try:
        return ill.loc[iid]['Name']
    except KeyError:
        #print ("Warning: %s not found in illumina sheet" % iid)
        return iid
        
def makeRow2(pos,row):
    meta = [row.CHR, pos, row.A2, row.A1, lookupRS(row['index']), 1, row.A1]
    freqs = [''] * (len(pops)-1) + [row.MAF]
    return sep.join(map(str, meta + freqs)) +'\n'

if len(sys.argv) == 1:
    for i in list(range(1,23)) + [25]: print(cmd % (i,i,i))
else:
    ## PLink input
    sep = '\t'
    bim = pd.read_csv("%s/AB_QC_2_2_3.bim"%wdir, sep='\t', index_col=1, header=None)
    frq = pd.read_csv("%s/AB_QC_2_2_3.frq"%wdir, sep='\s+', index_col='SNP')
    #frq file format, see: https://www.cog-genomics.org/plink2/formats#frq
    gme = pd.read_csv("%s/variome.trim_PanTro2_sampgenes.allsamples1.annot.tsv.gz" % gmeDir, sep='\t')
    #gme.columns = [colRename.get(col, col) for col in gme.columns.values]
    frq2 = pd.concat([bim, frq], axis=1) # merging based on SNP
    #np.histogram(frq["MAF"],  bins = 100) #plot this
    meta = 'CHROM POS REF ALT ID NALT MINOR'.split()
    pops = 'AFR AMR ASJ EAS FIN NFE OTH GME UAE'.split()

    ill = pd.read_csv("Data/infiniumomni5exome-4-v1-3-a1-b144-rsids.zip", sep='\t', index_col=0)

    ## TODO: figure out mapping between non-autosomal chromosome identifiers:
    ##       Ours
    ## Mt    23
    ## X     25
    ## Y     26
    chromosome = int(sys.argv[-1])
    chromosome1 = {25:'X'}.get(chromosome, chromosome)
    count = 0
    correctMatch = 0
    flippedMatch = 0
    f = open("Data/popAFs4_%s.csv"%chromosome, "w")
    f.write(sep.join(meta+pops)+"\n") 
    omniChr = frq2[frq2['CHR'] == chromosome].sort_values(3).reset_index().set_index(3)
    compatibilityFails = 0

    freqStats = {}
    ## GNOMAD annotation of our AFs
    # First starting with lines in Gnomad vcf, as this is a vcf flatfile, 
    # our genotype array is put into a dataframe

    gme22 = gme[gme['chrom'] == chromosome].set_index('pos') 
    gme22 = gme22[~gme22.index.duplicated(keep='first')] ## removing duplicate entries (wrt. index)
    drops = []
    for line in bgzf.open("%s/gnomad.genomes.r2.0.2.sites.chr%s.vcf.gz" % (gnomadDir, chromosome1)):
        line = line.decode() ## dealing with binary stuff
        if line.startswith("#"): continue
        gnomadLine = GnomADLine(line) 
        if gnomadLine.indel: continue
        whichAlt = 0 ## normally just consider first Alternative allele
        uae = False
        fs = None
        if gnomadLine.pos in omniChr.index:
            drops.append(gnomadLine.pos)
            fs = FreqStats(omniChr.loc[gnomadLine.pos])
            if fs.compatibilityCheck(gnomadLine):
                uae = True
                whichAlt = fs.whichAlt ## a bit hackish: whichAlt is wrt gnomad here

        gnomadLine.info2popfreqs(whichAlt)
        if uae:
            count += 1
            gnomadLine.popfreqs.append(('AF_UAE', fs.maf))

        if gnomadLine.pos in gme22.index:
            gmeAnnot = Annotation(gme22.loc[gnomadLine.pos])
            #fs.compatibilityCheck(gmeAnnot): ideally double check 
            if gmeAnnot.ref == gnomadLine.ref and gmeAnnot.alts[0] == gnomadLine.alts[0]:
                gnomadLine.popfreqs.append(('AF_GME', gmeAnnot.frq))
        f.write(makeRow(gnomadLine, fs))

    uaeNoGnomad = omniChr.drop(drops)
    print("UAE variants not in Gnomad", uaeNoGnomad.shape)
    for pos, row in uaeNoGnomad.iterrows():
        f.write(makeRow2(pos, row))
    f.close()

## a bit bizare results: ASJ is often the closest allele freq, followed by AMR ?!?!?! (ok no SouthAsian, but still) ->
## story of migration: Arabs ->  Andalucia ->