"vcftools --vcf uae_hgdp1LD.vcf --weir-fst-pop hgdp_North_Africa.txt --weir-fst-pop hgdp_Subsaharian_Africa.txt --weir-fst-pop hgdp_Central_South_Asia.txt --weir-fst-pop hgdp_Middle_Est.txt --out Fst_ME_SSA_NA_CSA"

import glob, subprocess
import pandas

hgdp glob.glob('../Data/hgdp_*.txt')
