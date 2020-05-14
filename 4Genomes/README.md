Various tools for performing genome analysis.

## Concordance

concordance.py: run without parameters to see help.
It requires files to be in plink (ped/map) format.
Concordance is calculated by requiring a number of criteria:
* Biallelic loci only
* variant present in dbSNP (SQL connection to dbSNP required)
* no strand confusion (if ped file has both major and minor allel being complementary to reference genome)
* no missingness

For more details see the doc string


## parseGnomad.py

Python script that goes through GnomAD in vcf format, extracts allele frequencies and adds them to the info column of the input vcf file. Performs a few sanity checks.

## table2vcf.py

convenience tool to transform data.

## snpsiftAnnotate.sh
Shell script that documents how snpSift was run

## vtools (aka variant tools)
Its usage is demonstrated in the dedicated README file.
