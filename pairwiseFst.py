import glob, sys
import subprocess
import os


pops = sorted(glob.glob('uaeIDs_*.txt')) + sorted(glob.glob('hgdp_*.txt'))
job = int(sys.argv[-1]) - 1
idx = [(k,l) for k in range(len(pops)) for l in range(k)]
#i,j = idx[job]
i,j = job,job
cmd = f'vcftools --vcf uae_hgdp1LD_fix.vcf --weir-fst-pop {pops[i]} --weir-fst-pop {pops[j]} --out pairFst_{i}_{j}'
print(cmd)
subprocess.call(cmd, shell=True)
        
