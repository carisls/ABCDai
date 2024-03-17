###############################################################
### Copyright©2024. Caris MPI, Inc. All rights reserved.###
###############################################################

#!/usr/bin/env python
# coding: utf-8


import os
import sys
import pysam
import subprocess
import yaml
from datetime import datetime
import shutil

import argparse
import pathlib
import traceback



launchCmd = " ".join(sys.argv)
app=sys.argv[0]
binDir=os.path.dirname(app)
if(binDir == ""):  binDir = "."

parser = argparse.ArgumentParser('Generate Mutationome features from BAM file using LoFreq and Sentieon TNHaplotyper2 for variant calling')
parser.add_argument("--bam", type=pathlib.Path, help="Input BAM filename and path", required=True)
parser.add_argument("-f", "--outputFolder", type=pathlib.Path, help="Output Folder", required=True)
parser.add_argument("-c", "--config", type=pathlib.Path, help="Configuration File", required=True)
parser.add_argument('--overwrite', action='store_true')
args = parser.parse_args()

sorted_bam = args.bam.as_posix()
outputFolder = args.outputFolder.as_posix()
configPath = args.config.as_posix()
overwrite = args.overwrite
naType = 'DNA'

if outputFolder[-1] == '/':
    outputFolder = outputFolder[:-1]

config = yaml.full_load(open(configPath,'rt'))
sentiEnv = config['sentiEnv']
sentiPath = config['sentiPath']
reference = config['hg38_RefSeq_file']
threads = config['threads']
pythonCmd = config['pythonCmd']
cprvUniverse = config['cprvUniverse']
binMeta = config['binMeta']
common_CH_FNAP = config['Common_CH']

GATK = config['gatk']
TARGETS_PADDED = config['targets_padded']
SNPEFF = config['snpeff']
JAVA = config['java']

os.makedirs(outputFolder, exist_ok=True)

def run_command(cmd, label):
    try:
        print(cmd)
        out = subprocess.check_output(cmd, shell=True)
        out = out.decode('ascii')
        print(out)
    except subprocess.CalledProcessError as e:
        print(label+" ERROR: ",e)
        print(e.output)

def call_LoFreq(bam, merged, overwrite):
    if not os.path.exists(merged) or (overwrite):
        fullvcf = merged.replace('.vcf','.full.vcf')
        cmd = f'{pythonCmd} scripts/lofreq_parallel.py -b {bam} -v {fullvcf} -c {configPath}'
        run_command(cmd, 'LoFreqParallel')
        
        cmd = f'{pythonCmd} scripts/VcfTrimmer.py --input {fullvcf} --output {merged} --reads 1 --force'
        run_command(cmd, 'LoFreqTrim_1_read')

        cmd = f'{GATK} IndexFeatureFile -I {merged}'
        run_command(cmd, 'LoFreqTrim_1_read_Index')


def call_mutect2(vcf, merged, sorted_bam, overwrite):
    if not os.path.exists(vcf) or (overwrite):
        sample_name = sorted_bam.split('/')[-1].split('.')[0]
        cmd = f'{sentiEnv}; time {sentiPath}/sentieon driver             -r {reference}             -t 24             --interval {TARGETS_PADDED}             -i {sorted_bam}             --algo TNhaplotyper2             --tumor_sample {sample_name}             --min_init_tumor_lod             -100             --min_tumor_lod -100             --callable_depth 0     --trim_soft_clip    --given {merged}     {vcf}'
        run_command(cmd, 'Mutect2')


def filter_mutect2(vcf, filtered_vcf, overwrite):
    if not os.path.exists(filtered_vcf) or (overwrite):
        fixedvcf = vcf.replace(".vcf",".fixed.vcf")
        fout = open(fixedvcf, "w")
        
        with open(binDir + "/header.vcf","r") as f:
            header=f.readlines()
           
        for line in header:
            fout.write(line)
            
        write=False
        with open(vcf,"r") as f:
            for line in f:
                if line.startswith("##SentieonCommandLine"):
                    write=True
                if write:
                    fout.write(line)
        fout.close()
        
        os.system("/bin/mv %s %s" % (vcf, vcf+".bak"))
        os.system("/bin/mv %s %s" % (fixedvcf, vcf))

        
        cmd = f'{GATK} FilterMutectCalls -V {vcf} -O {filtered_vcf} -R {reference}'
        run_command(cmd, 'FilterMutectCalls')

def call_left_align_trim_variants(filtered_vcf, split_vcf, overwrite):
    if not os.path.exists(split_vcf) or (overwrite):
        cmd = f'{GATK} LeftAlignAndTrimVariants             -R {reference}             -V {filtered_vcf}             -O {split_vcf}             --split-multi-allelics             --dont-trim-alleles             --keep-original-ac'
        run_command(cmd, 'LeftAlignAndTrimVariants')

def snpeff(split_vcf, anno_vcf, overwrite):
    if not os.path.exists(anno_vcf) or (overwrite):
        cmd_snpeff = f'{JAVA} -Xmx8g -jar {SNPEFF}/snpEff.jar         -q         -c {SNPEFF}/snpEff.config         -noStats hg38_20190May {split_vcf} > {anno_vcf}'
        run_command(cmd_snpeff, 'snpEff 4.3')

def filterVcf(anno_vcf, final_filtered_vcf, overwrite):
    if not os.path.exists(final_filtered_vcf) or (overwrite):
        cmd = f'{pythonCmd} scripts/FilterVcf.py --input {anno_vcf} --cprv {cprvUniverse} --output {final_filtered_vcf}'
        run_command(cmd, 'Filter_VCF')




bamFn = args.bam.name
sampleName = bamFn.split('.')[0].split('_S')[0]
merged = outputFolder+'/'+bamFn.replace('.bam', '.lofreq.vcf')
vcf = outputFolder+'/'+bamFn.replace('.bam', '.lofreq.unfiltered.vcf')
filtered_vcf = vcf.replace('unfiltered', 'filtered')
split_vcf = filtered_vcf.replace('.vcf','.split.vcf')
anno_vcf = split_vcf.replace('.vcf','.annotated.vcf')
final_filtered_vcf = anno_vcf.replace('.vcf','.filtered.vcf')
mutationome_out = f'{outputFolder}/{sampleName}.mutationome.csv'

call_LoFreq(sorted_bam, merged, overwrite)
call_mutect2(vcf, merged, sorted_bam, overwrite)
filter_mutect2(vcf, filtered_vcf, overwrite)
call_left_align_trim_variants(filtered_vcf, split_vcf, overwrite)
snpeff(split_vcf, anno_vcf, overwrite)
filterVcf(anno_vcf, final_filtered_vcf, overwrite)


print("Job Completed")





