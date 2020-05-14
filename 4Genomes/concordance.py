"""

usage: python concordance2.py <genotype.ped> <ngs.ped>

(assumes map files with same basename to be existent)
Concordance calculation between NGS and Omni5  Illumina Bead Chips
Expects two sets of map/ped files. The NGS version is derived from the vcf 
file with filter=PASS for all variants, using vcftools: 

vcftools --vcf 12742_recalibrated_variants_PASS.vcf --plink --out 12742_PASS

Needs MySQL connection to MySQL server hosted at BTS (localhost) in order to compare 
"""

from collections import Counter
import pdb
import sys, os

com= {'G':'C', 'C':'G', 'A':'T', 'T':'A'}

def complement(l1, l2):
    return (l1[0] == com[l2[0]] and l1[1] == com[l2[1]]) or (l1[1] == com[l2[0]] and l1[0] == com[l2[1]])
def disagreement(l1, l2):
    return sum((Counter(l1) - Counter(l2)).values())
def indel(g):
    return len(g[0]) > 1 or len(g[1]) > 1
def missing(g):
    return g[0] in '0' or g[1] == '0'
def undef(g):
    return not g[0] in ['A', 'C', 'G', 'T'] or not g[1] in ['A', 'C', 'G', 'T']
def concurWref(g, ref):
    return g[0] in [ref['ALT'], ref['REF']] and g[1] in [ref['ALT'], ref['REF']]
def totalagreement(reflookup, pos2gt_GT, pos2gt_NGS):
    """Main function to determine concordance between NGS sequencing (vcf->ped) and genotype array data from
     the same subject. To be as precise as possible, a number of cases are ruled out for comparison: 
    1. indels (not called on chip)
    2. missingness in either dataset
    3. if beadchip alleles are not in accordance with dbSNP
    4. potential strandconfusion, e.g. genotype array alleles are the exact complement 
    5. Multiallelic loci
"""
    disagreements = 0
    agreements = 0
    missingcount = 0
    indelcount = 0
    strandconfusion = 0
    beadChipErrors = 0
    multiallelic = 0
    discordance = 0
    
    for (chrom, pos), referenceList in reflookup.items():
        reference = set(referenceList)
        ##check for possible strand confusion
        g1 = pos2gt_GT[(chrom, str(pos))]
        if not (chrom, str(pos)) in pos2gt_NGS:
            print("This should not happen")
            continue
            #g2 = [referenceList[0], referenceList[0]] 
        else:
            g2 = pos2gt_NGS[(chrom, str(pos))]
        if indel(g2): 
             indelcount += 1
             continue
        if missing(g1) or missing(g2): 
             missingcount += 1
             continue
        
        if not reference.issuperset(g1): 
            beadChipErrors += 1
            continue

        da = disagreement(g1, g2)
        if da != 0:
            excuse = False
            if complement(g1, g2):
                strandconfusion += 1
                excuse = True
            if len(reference) > 2: 
                multiallelic += 1
                excuse = True
            if not excuse:
                discordance += 1
                #print ('\t'.join(g1 + g2 + [chrom, str(pos), str(da)] + referenceList))
                # ['G', 'G'] ['G', 'A'] ('3', 16409491) 1
            disagreements += 1
        else:
            agreements += 1

    ##Final report
    print ("Comparison based on:                            %s SNPs"% (len(reflookup)))
    print ("Total disagreements:                            %s (incl. pot. strand confusion, Multi-Allelic Loci)"% disagreements)
    print ("Missing:                                        %s" % missingcount)
    print ("Indels:                                         %s" % indelcount)
    print ("Bead Chip alleles not a subset of reference:    %s" % beadChipErrors)
    print ("Agreements:                                     %s" % agreements)
    print ("Strandconfusion:                                %s" % strandconfusion)
    print ("Multiallelic loci:                              %s" % multiallelic)
    print ("Discordance:                                    %s" % discordance)
    print ("Concordance not ruling out strandconfusion etc. %.4f%%" % 
           (100.*agreements/(disagreements + agreements)))
    print ("Concordance (cleaned)                           %.4f%%" % (100.*agreements/(discordance + agreements)))

# see file formats http://www.gwaspi.org/?page_id=145
def parsePED(pedfile, mapfile, sep='\t'):
    '''simultaneously parsing ped and map file'''
    def mapping(line):
        fields = line.strip().split()
        return (fields[0], fields[-1])

    def readcolumns(tsv, maxtabs=2):
        tabcount = 0
        info = ''
        while tabcount < maxtabs:
            c = tsv.read(1)
            info += c
            if not c: break
            if re.match('\s', c): tabcount += 1
        return re.split('\s+', info)[:maxtabs]

    pedf = open(pedfile)
    famID, sampleID = readcolumns(pedf, maxtabs=6)[:2]
    ## creating pos2gt = {('21', '84932002'): ['C', 'G'], ...}
    pos2gt = dict([( (line.strip().split()[0], line.strip().split()[-1]) , readcolumns(pedf)) for line in open(mapfile)])
    unparsedBytes = len(pedf.read()) # reads the rest, if any (shouldn't be)
    pedf.close()
    return pos2gt

def singletons(variants):
    return [variant for variant in variants.split(',') if len(variant)==1]

def flatten(lol):
    return [y for x in lol for y in x]
def referenceForOverlap(overlap):
    import copy
    lookup = {}
    for (chrom, pos) in overlap:
        cursor.execute('SELECT REF, ALT, INFO, GMAF FROM annotator.dbSNP WHERE CHR="%s" and POS=%s AND INFO="SNV"' % (chrom, pos))
        data = cursor.fetchall()
        ref = list(set(flatten([singletons(rec['REF']) for rec in data])))
        if len(ref) > 1: print ("Warning, multiple REF entry:", (chrom, pos))
        alt = flatten([singletons(rec['ALT']) for rec in data])
        if ref or alt: lookup[(chrom, pos)] = ref + alt
        #for rec in data:
        #    lookup[(chrom, pos)] = copy.copy(data[0])
        #    
        #    if len(data) > 1:
        #        print "Multiple results for Chr%s Pos %s" % (chrom, pos),
        #        print data[:3]
    return lookup

if __name__ == '__main__':
    import mysql.connector
    import re
    import sys

    conn = mysql.connector.connect(user='root', password="P@ssw0rdku123", database="annotator", host="localhost")
    cursor = conn.cursor(dictionary=True)

    if len(sys.argv) < 2:
         print(__doc__)
         sys.exit(1)
         # was 12742Genotype.ped, the full unfiltered (no QC) file
         #gtpedfile = '/home/ahenschel/Concordance/Data/mergedQCLD_10745.ped'  
         #ngspedfile = '/bmshare/gihan/RawData/10745/10745_PASS_genotype.ped'
    else:
         gtpedfile = sys.argv[1]
         gtmapfile = os.path.splitext(gtpedfile)[0] + '.map'
         ngspedfile = sys.argv[2]
         ngsmapfile = os.path.splitext(ngspedfile)[0] + '.map'


    geno = [(line.split()[0], int(line.split()[-1])) for line in open(gtmapfile)]
    ngs = [(line.split()[0], int(line.split()[-1])) for line in open(ngsmapfile)]
    print ("parsing map files done")
    overlap = set(geno) & set(ngs)
    try:
        reflookup = referenceForOverlap(overlap)
        print (len(reflookup))
    finally:
        conn.close()

    ## get all (chromosome, position) tuples for which we have data from both technologies
    # 252037 for QC dataset
    pos2gt_GT  = parsePED(gtpedfile, gtmapfile, sep=' ')
    pos2gt_NGS = parsePED(ngspedfile, ngsmapfile)

    sc = totalagreement(reflookup, pos2gt_GT, pos2gt_NGS)
'''
Steps to calculate concordance
plink --bfile CVRL --make-bed --maf 0.01 --geno 0.01 --hwe 0.001 --out CVRL_QC_2_2_3
cat > keep10187.txt
59 9239299151_R03C01
plink --bfile CVRL_QC_2_2_3 --keep keep10187.txt --recode tab --out CVRL_10187

## NGS
/bmshare/gihan/results/10187> vcftools --vcf 10187_snps_passed.vcf --plink --out ~/Concordance/Data/10187_PASS

run concordance2.py Data/CVRL_10187.ped Data/10187_PASS.ped
Comparison based on: 332968 SNPs
Total disagreements: 5271
Agreements:          167025
Strandconfusion:     5130
Multiallelic loci:   4886
Discordance:         111
Concordance not ruling out strandconfusion etc.  96.9407%
Concordance 99.9336%
'''
