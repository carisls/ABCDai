###############################################################
### Copyright©2024. Caris MPI, Inc. All rights reserved.###
###############################################################


import os
import sys
import yaml
import subprocess

from joblib import Parallel, delayed

import argparse
import pathlib

# In[7]:

parser = argparse.ArgumentParser('Runs LoFreq in 12 parallel threads using 12 bed files to cover the entire genome.')
parser.add_argument('-b',"--bam", type=pathlib.Path, help="Input BAM filename and path", required=True)
parser.add_argument('-v',"--vcf", type=pathlib.Path, help="Output VCF filename and path", required=True)
parser.add_argument("-c", "--config", type=pathlib.Path, help="Configuration File", required=True)
args = parser.parse_args()

inputBam = args.bam.as_posix()
outputVCF = args.vcf.as_posix()
configPath = args.config.as_posix()


config = yaml.full_load(open(configPath))

hg38_fa = config['hg38_RefSeq_file']
seq_dict = config['seq_dict']
GATK = config['gatk']

def run_command(cmd):
    try:
        print(cmd)
        out = subprocess.check_output(cmd, shell=True)
        out = out.decode('ascii')
        print(out)
    except subprocess.CalledProcessError as e:
        print("ERROR: ",e)
        print(e.output)   


cmds = []
vcfShards = []
for i in range(1,13):
    vcfShard = outputVCF+f'.GenomeShard_{i}Of12.vcf'
    cmd  = f"lofreq call -f {hg38_fa} "
    cmd += f"-l data/GenomeShard_{i}Of12.bed --no-default-filter --no-baq --no-mq --sig 1 --call-indels "
    cmd += f"--bonf 1  -o {vcfShard} {inputBam}"
    cmds.append(cmd)
    vcfShards.append(vcfShard)
    
    
Parallel(n_jobs = len(cmds))(delayed(run_command)(cmd) for cmd in cmds)

inputVcfs=''
for vcf in vcfShards:
    inputVcfs += f' I={vcf}'

cmd = f'{GATK} MergeVcfs {inputVcfs} O={outputVCF} D={seq_dict}'
run_command(cmd)

for vcf in vcfShards:
    os.remove(vcf)

print("Job Completed")


# In[ ]:




