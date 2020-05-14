
## vtools usage
## create a project popAFs
vtools init popAFs #-> popAFs.proj
## Make sure there is a popAF.fmt file so to tell vtools about the format
 2041  cp /research/gutsybugs/KUMI/PopAFs/variant_tools/format/* .
 2043  cd /research/gutsybugs/KUMI/PopAFs
## After project creation, import vcf files
 2053  vtools execute import_vcf --input popAFs*vcf.gz --output popAF.fmt --build hg19
 2054  vtools show tables

## Make tables
 2061  vtools select variant "AF_Z>4" -t highZ4
 2067  vtools select variant "AF_Z<-4" -t lowZ4

## Run snpEff, only on highZ4 and lowZ4 tables
 2062  vtools execute snpEff eff --snpeff_path /home/shared/snpeff/snpEff --var_table highZ4
 2069  vtools execute snpEff eff --snpeff_path /home/shared/snpeff/snpEff --var_table lowZ4

## vtool commands that create tabular output such as 10745_Z-1.96_clnvarNotBenign.csv in Data/4genomes/Results
#grep "select variant" popAF.log |grep "clinvar.CLNSIG" | grep AF_UAE

select variant AF_Z>4 "clinvar.CLNSIG LIKE '%4%' OR clinvar.CLNSIG LIKE '%5%' OR clinvar.CLNSIG LIKE '%0%'" -o variant_id chr pos ref alt variant.name AF_Z AF_UAE AF_GME AF_gMAX clinvar.CLNSIG clinvar.CLNDBN clinvar.CLNDSDB

select variant AF_Z>4 "clinvar.CLNSIG LIKE '%4%' OR clinvar.CLNSIG LIKE '%5%' OR clinvar.CLNSIG LIKE '%0%'" -o variant_id chr pos ref alt variant.name AF_Z AF_UAE AF_GME AF_gMAX clinvar.CLNSIG clinvar.CLNDBN clinvar.CLNDSDB --header

select variant AF_Z>4 "clinvar.CLNSIG LIKE '%4%' OR clinvar.CLNSIG LIKE '%5%' OR clinvar.CLNSIG LIKE '%0%'" -o variant_id chr pos ref alt variant.name AF_Z AF_UAE AF_GME AF_AMR AF_EAS AF_ASJ AF_AFR AF_NFE AF_gMAX clinvar.CLNSIG clinvar.CLNDBN clinvar.CLNDSDB --header

select variant AF_Z<-4 "clinvar.CLNSIG LIKE '%4%' OR clinvar.CLNSIG LIKE '%5%' OR clinvar.CLNSIG LIKE '%0%'" -o variant_id chr pos ref alt variant.name AF_Z AF_UAE AF_GME AF_AMR AF_EAS AF_ASJ AF_AFR AF_NFE AF_gMAX clinvar.CLNSIG clinvar.CLNDBN clinvar.CLNDSDB --header

select variant AF_Z>2 "clinvar.CLNSIG LIKE '%4%' OR clinvar.CLNSIG LIKE '%5%' OR clinvar.CLNSIG LIKE '%0%'" -o variant_id chr pos ref alt variant.name AF_Z AF_UAE AF_GME AF_AMR AF_EAS AF_ASJ AF_AFR AF_NFE AF_gMAX clinvar.CLNSIG clinvar.CLNDBN clinvar.CLNDSDB --header

2019-02-14 17:01:46,725: DEBUG: select variant AF_Z>1.96 "clinvar.CLNSIG LIKE '%4%' OR clinvar.CLNSIG LIKE '%5%' OR clinvar.CLNSIG LIKE '%0%'" -o variant_id chr pos ref alt variant.name AF_Z AF_UAE AF_GME AF_AMR AF_EAS AF_ASJ AF_AFR AF_NFE AF_gMAX clinvar.CLNSIG clinvar.CLNDBN clinvar.CLNDSDB --header

## in ipython (also calling vtools):
##Create result stats about snpeff annotations of loci with sign. diff. af's
import pandas as pd

!vtools select highZ4 -o chr pos ref variant_id alt AF_UAE AF_gMIN AF_Z EFF > highZ4_EFF1.csv
df=pd.read_csv('highZ4_EFF1.csv', delimiter='\t', header=None)

!vtools select highZ4 "EFF LIKE 'STOP_GAINED%'" -o chr pos ref variant_id alt AF_UAE AF_gMIN AF_Z EFF -l5
!vtools select highZ4 "EFF LIKE 'START_GAINED%' OR EFF LIKE 'STOP_GAINED%'" -o chr pos ref variant_id alt AF_UAE AF_gMIN AF_Z EFF -l5
!vtools select highZ4 "EFF LIKE 'START_GAINED%' OR EFF LIKE 'STOP_GAINED%' OR EFF LIKE 'NON_SYNONYMOUS_CODING%'" -o chr pos ref variant_id alt AF_UAE AF_gMIN AF_Z EFF -t StartStopNonsyn
!vtools select highZ4 "EFF LIKE 'START_GAINED%' OR EFF LIKE 'STOP_GAINED%' OR EFF LIKE 'NON_SYNONYMOUS_CODING%'" -t StartStopNonsyn
!vtools select StartStopNonsyn 'EFF LIKE 'START_GAINED%' -o chr pos ref variant_id alt AF_UAE AF_gMIN AF_Z EFF
!vtools select highZ4 'EFF LIKE 'START_GAINED%' -o chr pos ref variant_id alt AF_UAE AF_gMIN AF_Z EFF
!vtools select StartStopNonsyn "EFF LIKE 'START_GAINED%'" -o chr pos ref variant_id alt AF_UAE AF_gMIN AF_Z EFF
!vtools select StartStopNonsyn --count
!vtools show table variant


!vtools select StartStopNonsyn -o variant_id bin chr pos ref alt name AF_AFR AF_AMR AF_ASJ AF_EAS AF_FIN AF_NFE AF_OTH AF_GME AF_UAE AF_gMAX AF_gMIN AF_Z AF_Z1 EFF -l5
!vtools select StartStopNonsyn -o variant_id bin chr pos ref alt name AF_AFR AF_AMR AF_ASJ AF_EAS AF_FIN AF_NFE AF_OTH AF_GME AF_UAE AF_gMAX AF_gMIN AF_Z AF_Z1 EFF --header > highZ4_EFFdelet.csv
!vtools select lowZ4 "EFF LIKE 'START_GAINED%' OR EFF LIKE 'STOP_GAINED%' OR EFF LIKE 'NON_SYNONYMOUS_CODING%'" -t StartStopNonsynLow
!vtools select lowZ4 "EFF LIKE 'START_GAINED%' OR EFF LIKE 'STOP_GAINED%' OR EFF LIKE 'NON_SYNONYMOUS_CODING%'" -t StartStopNonsynLow
!vtools select StartStopNonsynLow -o variant_id bin chr pos ref alt name AF_AFR AF_AMR AF_ASJ AF_EAS AF_FIN AF_NFE AF_OTH AF_GME AF_UAE AF_gMAX AF_gMIN AF_Z AF_Z1 EFF --header >> highZ4_EFFdelet.csv

df=pd.read_csv('highZ4_EFFdelet.csv', delimiter='\t', header=None)


## Further selection using pandas algebra
