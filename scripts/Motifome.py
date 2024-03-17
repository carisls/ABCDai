###############################################################
### Copyright©2024. Caris MPI, Inc. All rights reserved.###
###############################################################

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

parser = argparse.ArgumentParser('Generate Motifome features from BAM file using LoFreq and Sentieon TNHaplotyper2 for variant calling')
parser.add_argument("--bam", type=pathlib.Path, help="Input BAM filename and path", required=True)
parser.add_argument("-f", "--outputFolder", type=pathlib.Path, help="Output Folder", required=True)
parser.add_argument("-c", "--config", type=pathlib.Path, help="Configuration File", required=True)
parser.add_argument('--overwrite', action='store_true')
args = parser.parse_args()

sorted_bam = args.bam.as_posix()
outputFolder = args.outputFolder.as_posix()
configPath = args.config.as_posix()
overwrite = args.overwrite

if outputFolder[-1] == '/':
    outputFolder = outputFolder[:-1]

config = yaml.full_load(open(configPath,'rt'))
pythonCmd = config['pythonCmd']
inclusionRegionBed = config['cbr_bed']
exclusionRegionBed = config['contamination_bed']

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

        
        
        
bamFn = args.bam.name
sampleName = bamFn.split('.')[0].split('_S')[0]
if not os.path.exists(f'{outputFolder}/{bamFn}.FragLen.pkl.gz') or (overwrite):
    cmd = f'{pythonCmd} scripts/get_fragment_metrics_singlethread_blacklist.py -i {sorted_bam} -o {outputFolder} -c {inclusionRegionBed} -b {exclusionRegionBed}; {pythonCmd} scripts/motif_feature_extraction.py -i {outputFolder}/{bamFn}.FragLen.pkl.gz -o {outputFolder}'
    run_command(cmd, 'Motifome')    
    
print("Job Completed")
