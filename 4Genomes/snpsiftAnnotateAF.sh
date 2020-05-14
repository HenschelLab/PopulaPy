## Parallelizing shell script to annotate allele frequency vcf files and previous (main) vcf files  

suffix="_snpeff_annotated_gatk_clinvar2_dbsnp2017071_vartype_dbNSFP"
snpsift="/usr/share/BTSSRV/snpEff/SnpSift.jar"
targetDir="/bmshare/ahenschel/AnnotatedGenomes/$1/"

task() {
    sleep $(($2 + ($3-1)*60));
    #echo "$1 $2 $3"
    nohup java -Xmx10g -jar ${snpsift} annotate /home/ahenschel/KUMI/Data/popAFs3_$3.vcf.gz /bmshare/gihan/RawData/$1/Chroms/$1${suffix}_chr$3a.vcf.gz > ${targetDir}/$1${suffix}_chr${3}.
vcf 2> Nohup/snpSift_$1_$3.nohup.out
}

for i in {1..22}; do 
    task $1 $2 $i &
done